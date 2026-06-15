/*****************************************************************************
 *   Copyright (c) 2020, Hobu, Inc. (info@hobu.co)                           *
 *   Modified by Kyle Kam (kylekam@mach9.io)                                 *
 *                                                                           *
 *   All rights reserved.                                                    *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 ****************************************************************************/

#include <algorithm>
#include <numeric>

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

#include <pdal/PDALUtils.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/filters/SortFilter.hpp>

#include <lazperf/lazperf.hpp>
#include <lazperf/writers.hpp>

#include "ChunkWriter.hpp"
#include "PyramidManager.hpp"

#include <stringconv.hpp>  // untwine/os

namespace untwine
{
namespace bu
{

ChunkWriter::ChunkWriter(PyramidManager& manager, const BaseInfo& b) :
    m_manager(manager), m_b(b), m_pool(NumChunkWriters, ChunkWriterQueueSize)
{
    m_pool.trap(true, "Unknown error in ChunkWriter");

#ifndef _WIN32
    // Late materialization: the attribute store was fully written during EPF
    // (BU is constructed after EPF completes).
    if (m_b.opts.lateMaterialization)
    {
        std::string attrPath = m_b.opts.tempDir + "/" + AttributeStoreFilename;
        m_attrFd = ::open(attrPath.c_str(), O_RDONLY);
        if (m_attrFd < 0)
            throw FatalError("Can't open attribute store '" + attrPath + "'.");
    }
#endif
}

ChunkWriter::~ChunkWriter()
{
    closeAttrStore();
}

void ChunkWriter::closeAttrStore()
{
#ifndef _WIN32
    if (m_attrFd >= 0)
    {
        ::close(m_attrFd);
        m_attrFd = -1;
    }
#endif
}

void ChunkWriter::enqueue(Chunk&& chunk)
{
    // The task must be copyable to live in a std::function, so hold the chunk
    // through a shared pointer.
    auto c = std::make_shared<Chunk>(std::move(chunk));
    m_pool.add([this, c]()
    {
        write(*c);
    });
}

void ChunkWriter::stop()
{
    // join() drains the queue - workers only exit once all queued chunks have
    // been written.
    m_pool.join();
    closeAttrStore();
    StringList errors = m_pool.clearErrors();
    if (errors.size())
        throw FatalError(errors.front());
}

// Late materialization: the view arrives holding only X/Y/Z; fill in every other
// dimension from the attribute store. Must run before anything reads non-xyz dims:
// stats (write()), the GpsTime sort and compression (createChunk()).
void ChunkWriter::gather(Chunk& chunk)
{
#ifndef _WIN32
    if (m_attrFd < 0 || chunk.view->size() == 0)
        return;

    const uint64_t W = (uint64_t)m_b.attrRecordSize;

    // Visit view points in id order so reads of the input-ordered attribute store
    // coalesce into sequential runs (a voxel's ids cluster - scan data is coherent).
    std::vector<uint32_t> order(chunk.ids.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
        [&ids = chunk.ids](uint32_t a, uint32_t b){ return ids[a] < ids[b]; });

    // Resolve the non-xyz dims against this chunk's table layout. X/Y/Z occupy
    // wide-record bytes [0, SlimIdOffset); every other dim's attribute-store offset
    // is its wide offset minus that prefix.
    struct GatherDim
    {
        pdal::Dimension::Id dim;
        pdal::Dimension::Type type;
        int offset;
    };
    std::vector<GatherDim> gdims;
    for (const FileDimInfo& fdi : m_b.dimInfo)
        if (fdi.offset >= SlimIdOffset)
            gdims.push_back({chunk.view->layout()->findDim(fdi.name), fdi.type,
                fdi.offset - SlimIdOffset});

    std::vector<char> buf;
    size_t i = 0;
    while (i < order.size())
    {
        // Coalesce consecutive ids into one read.
        size_t j = i + 1;
        while (j < order.size() && chunk.ids[order[j]] == chunk.ids[order[j - 1]] + 1)
            j++;

        uint64_t firstId = chunk.ids[order[i]];
        size_t want = (j - i) * W;
        buf.resize(want);
        size_t got = 0;
        while (got < want)
        {
            ssize_t r = ::pread(m_attrFd, buf.data() + got, want - got,
                (off_t)(firstId * W + got));
            if (r <= 0)
                throw FatalError("Failure reading attribute store.");
            got += r;
        }

        for (size_t k = i; k < j; ++k)
        {
            char *src = buf.data() + (k - i) * W;
            for (const GatherDim& gd : gdims)
                chunk.view->setField(gd.dim, gd.type, order[k], src + gd.offset);
        }
        i = j;
    }
#endif
}

void ChunkWriter::write(Chunk& chunk)
{
    gather(chunk);

    // For single file output we need the stats for
    if (m_b.opts.stats)
    {
        for (pdal::PointId id = 0; id < chunk.view->size(); ++id)
        {
            for (auto& sp : chunk.stats)
            {
                pdal::Dimension::Id dim = sp.first;
                Stats& s = sp.second;
                s.insert(chunk.view->getFieldAs<double>(dim, id));
            }
        }
    }

    createChunk(chunk);
    m_manager.logOctant(chunk.key, (int)chunk.count, chunk.stats);
}

void ChunkWriter::sortChunk(pdal::PointViewPtr view)
{
    pdal::BufferReader r;
    r.addView(view);

    pdal::SortFilter s;
    s.setInput(r);
    pdal::Options o;
    o.add("dimension", "GpsTime");
    s.setOptions(o);

    s.prepare(view->table());
    s.execute(view->table());
}

void ChunkWriter::createChunk(const Chunk& c)
{
    using namespace pdal;

    if (c.view->size() == 0)
    {
        m_manager.newChunk(c.key, 0, 0);
        return;
    }

    // Sort the chunk on GPS time.
    if (c.view->layout()->hasDim(Dimension::Id::GpsTime))
        sortChunk(c.view);

    PointLayoutPtr layout = c.view->layout();

    int ebCount {0};
    for (DimType dim : c.extraDims)
        ebCount += layout->dimSize(dim.m_id);

    std::vector<char> buf(lazperf::baseCount(m_b.pointFormatId) + ebCount);
    lazperf::writer::chunk_compressor compressor(m_b.pointFormatId, ebCount);
    for (PointId idx = 0; idx < c.view->size(); ++idx)
    {
        PointRef point(*c.view, idx);
        fillPointBuf(point, buf, layout->findDim(UntwineBitsDimName), c.extraDims);
        compressor.compress(buf.data());
    }
    std::vector<unsigned char> chunk = compressor.done();

    uint64_t location = m_manager.newChunk(c.key, chunk.size(), (uint32_t)c.view->size());

    std::ofstream out(os::toNative(m_b.opts.outputName),
        std::ios::out | std::ios::in | std::ios::binary);
    out.seekp(std::ofstream::pos_type(location));
    out.write(reinterpret_cast<const char *>(chunk.data()), chunk.size());
    out.close();
    if (!out)
        throw FatalError("Failure writing to file '" + m_b.opts.outputName + "'.");
}

void ChunkWriter::fillPointBuf(pdal::PointRef& point, std::vector<char>& buf,
    pdal::Dimension::Id bitsDim, const pdal::DimTypeList& extraDims)
{
    using namespace pdal;

    LeInserter ostream(buf.data(), buf.size());

    bool hasTime = true; //  m_lasHeader.hasTime();
    bool hasColor = m_b.pointFormatId == 7 || m_b.pointFormatId == 8;
    bool hasInfrared = m_b.pointFormatId == 8;

    // we always write the base fields
    using namespace Dimension;

    uint8_t returnNumber(1);
    uint8_t numberOfReturns(1);
    if (point.hasDim(Id::ReturnNumber))
        returnNumber = point.getFieldAs<uint8_t>(Id::ReturnNumber);
    if (point.hasDim(Id::NumberOfReturns))
        numberOfReturns = point.getFieldAs<uint8_t>(Id::NumberOfReturns);

    auto converter = [](double d, Dimension::Id dim) -> int32_t
    {
        int32_t i(0);

        if (!Utils::numericCast(d, i))
            throw FatalError("Unable to convert scaled value (" +
                Utils::toString(d) + ") to "
                "int32 for dimension '" + Dimension::name(dim) +
                "' when writing LAS/LAZ file.");
        return i;
    };

    double x = point.getFieldAs<double>(Id::X);
    int32_t xi = converter((x - m_b.xform.offset.x) / m_b.xform.scale.x, Id::X);
    double y = point.getFieldAs<double>(Id::Y);
    int32_t yi = converter((y - m_b.xform.offset.y) / m_b.xform.scale.y, Id::Y);
    double z = point.getFieldAs<double>(Id::Z);
    int32_t zi = converter((z - m_b.xform.offset.z) / m_b.xform.scale.z, Id::Z);

    ostream << xi << yi << zi;

    ostream << point.getFieldAs<uint16_t>(Id::Intensity);
    ostream << (uint8_t)(returnNumber | (numberOfReturns << 4));
    ostream << point.getFieldAs<uint8_t>(bitsDim);
    ostream << point.getFieldAs<uint8_t>(Id::Classification);

    uint8_t userData = point.getFieldAs<uint8_t>(Id::UserData);
     // Guaranteed to fit if scan angle rank isn't wonky.
    int16_t scanAngleRank =
        static_cast<int16_t>(std::round(
            point.getFieldAs<float>(Id::ScanAngleRank) / .006f));
    ostream << userData << scanAngleRank;

    ostream << point.getFieldAs<uint16_t>(Id::PointSourceId);

    if (hasTime)
        ostream << point.getFieldAs<double>(Id::GpsTime);

    if (hasColor)
    {
        ostream << point.getFieldAs<uint16_t>(Id::Red);
        ostream << point.getFieldAs<uint16_t>(Id::Green);
        ostream << point.getFieldAs<uint16_t>(Id::Blue);
    }

    if (hasInfrared)
        ostream << point.getFieldAs<uint16_t>(Id::Infrared);

    Everything e;
    for (auto& dim : extraDims)
    {
        point.getField((char *)&e, dim.m_id, dim.m_type);
        Utils::insertDim(ostream, dim.m_type, e);
    }
}

} // namespace bu
} // namespace untwine
