/*****************************************************************************
 *   Copyright (c) 2020, Hobu, Inc. (info@hobu.co)                           *
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
#include <cmath>
#include <cstring>
#include <fstream>
#include <numeric>

#include <pdal/PDALUtils.hpp>
#include <pdal/util/Inserter.hpp>

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
    m_manager(manager), m_b(b),
    // Pool size and queue depth are env-tunable (defaults NumChunkWriters / ChunkWriterQueueSize)
    // so the compression thread budget can be matched to the box's I/O regime without a rebuild:
    // on RAM-constrained / slow-I/O boxes the build threads block on faults, so more writers fill
    // idle cores; on fast-storage boxes the build threads already saturate cores, so fewer writers
    // avoid oversubscribing. See untwine::envOverride.
    m_pool(envOverride("UNTWINE_NUM_CHUNK_WRITERS", NumChunkWriters),
           envOverride("UNTWINE_CHUNK_QUEUE_SIZE", ChunkWriterQueueSize)),
    // Build the packed-record layout straight from m_b.dimInfo, whose canonical ids/offsets are set
    // in prep/FilePrep::fillMetadata. Every packed-byte read below goes through the PackedPoint
    // views this produces; the ctor validates the layout once (fields in-bounds, non-overlapping).
    m_layout(m_b.dimInfo, m_b.pointSize)
{
    m_pool.trap(true, "Unknown error in ChunkWriter");

    resolveStdLocs();
}

// Resolve the standard fillPointBuf field locations from m_layout. Run-invariant (depends only on
// the layout, fixed for the whole job), so done once here rather than per chunk. Most resolve by
// canonical id; the packed-bits field has a dynamically-assigned id, so it's matched by name.
void ChunkWriter::resolveStdLocs()
{
    using namespace pdal;

    m_stdLocs.x = m_layout.field(Dimension::Id::X);
    m_stdLocs.y = m_layout.field(Dimension::Id::Y);
    m_stdLocs.z = m_layout.field(Dimension::Id::Z);
    m_stdLocs.intensity = m_layout.field(Dimension::Id::Intensity);
    m_stdLocs.returnNumber = m_layout.field(Dimension::Id::ReturnNumber);
    m_stdLocs.numberOfReturns = m_layout.field(Dimension::Id::NumberOfReturns);
    m_stdLocs.classification = m_layout.field(Dimension::Id::Classification);
    m_stdLocs.userData = m_layout.field(Dimension::Id::UserData);
    m_stdLocs.scanAngleRank = m_layout.field(Dimension::Id::ScanAngleRank);
    m_stdLocs.pointSourceId = m_layout.field(Dimension::Id::PointSourceId);
    m_stdLocs.gpsTime = m_layout.field(Dimension::Id::GpsTime);
    m_stdLocs.red = m_layout.field(Dimension::Id::Red);
    m_stdLocs.green = m_layout.field(Dimension::Id::Green);
    m_stdLocs.blue = m_layout.field(Dimension::Id::Blue);
    m_stdLocs.infrared = m_layout.field(Dimension::Id::Infrared);
    m_stdLocs.bits = m_layout.field(UntwineBitsDimName);
}

void ChunkWriter::enqueue(Chunk&& chunk)
{
    // The task must be copyable to live in a std::function, so hold the chunk through a shared
    // pointer.
    auto c = std::make_shared<Chunk>(std::move(chunk));
    m_pool.add([this, c]()
    {
        write(*c);
    });
}

void ChunkWriter::stop()
{
    // join() drains the queue; workers only exit once all queued chunks have been written.
    m_pool.join();
    StringList errors = m_pool.clearErrors();
    if (errors.size())
        throw FatalError(errors.front());
}

void ChunkWriter::write(Chunk& chunk)
{
    // Stats are accumulated here, on the worker thread, rather than on the build thread. Each
    // (id, Stats) pair keys on the same id the producer registered, so we resolve that id's Field
    // once and read the value straight from each packed record.
    if (m_b.opts.stats)
    {
        const char* base = chunk.records.data();
        for (auto& sp : chunk.stats)
        {
            pdal::Dimension::Id dim = sp.first;
            Stats& s = sp.second;
            Field loc = m_layout.field(dim);
            if (!loc.present())
                continue;
            for (size_t i = 0; i < chunk.count; ++i)
                s.insert(m_layout.point(base, i).get<double>(loc));
        }
    }

    createChunk(chunk);
    m_manager.logOctant(chunk.key, (int)chunk.count, chunk.stats);
}

void ChunkWriter::createChunk(const Chunk& c)
{
    using namespace pdal;

    if (c.count == 0)
    {
        m_manager.newChunk(c.key, 0, 0);
        return;
    }

    const char* base = c.records.data();

    // Extra-byte size from the dim types directly; no PointLayout needed. Each extra dim's Field
    // (offset + stored type) comes from the layout; the LAZ extra-byte layout uses dim.m_type,
    // which matches the stored type.
    int ebCount {0};
    std::vector<Field> extraLocs;
    extraLocs.reserve(c.extraDims.size());
    for (const DimType& dim : c.extraDims)
    {
        ebCount += (int)Dimension::size(dim.m_type);
        // Offset from the layout; type pinned to dim.m_type (the LAZ extra-byte type, == the stored
        // type), so an absent extra dim still emits Dimension::size(dim.m_type) zero bytes.
        Field eloc = m_layout.field(dim.m_id);
        eloc.type = dim.m_type;
        extraLocs.push_back(eloc);
    }

    // Sort on GPS time, replacing the old PDAL SortFilter. We compute a stable ascending order of
    // point indices by the GpsTime double read at its packed offset, then emit in that order. If
    // there's no GpsTime dimension, emit in natural order.
    std::vector<uint32_t> order(c.count);
    std::iota(order.begin(), order.end(), 0u);
    if (m_stdLocs.gpsTime.present())
    {
        const Field gps = m_stdLocs.gpsTime;
        const RecordLayout& layout = m_layout;
        std::stable_sort(order.begin(), order.end(),
            [base, gps, &layout](uint32_t a, uint32_t b)
            {
                return layout.point(base, a).get<double>(gps) <
                    layout.point(base, b).get<double>(gps);
            });
    }

    std::vector<char> buf(lazperf::baseCount(m_b.pointFormatId) + ebCount);
    lazperf::writer::chunk_compressor compressor(m_b.pointFormatId, ebCount);
    for (uint32_t idx : order)
    {
        fillPointBuf(m_layout.point(base, idx), buf, m_stdLocs, extraLocs);
        compressor.compress(buf.data());
    }
    std::vector<unsigned char> chunk = compressor.done();

    uint64_t location = m_manager.newChunk(c.key, chunk.size(), (uint32_t)c.count);

    std::ofstream out(os::toNative(m_b.opts.outputName),
        std::ios::out | std::ios::in | std::ios::binary);
    out.seekp(std::ofstream::pos_type(location));
    out.write(reinterpret_cast<const char *>(chunk.data()), chunk.size());
    out.close();
    if (!out)
        throw FatalError("Failure writing to file '" + m_b.opts.outputName + "'.");
}

void ChunkWriter::fillPointBuf(const PackedPoint& pt, std::vector<char>& buf, const StdLocs& loc,
    const std::vector<Field>& extraLocs)
{
    using namespace pdal;

    LeInserter ostream(buf.data(), buf.size());

    bool hasTime = true;
    bool hasColor = m_b.pointFormatId == 7 || m_b.pointFormatId == 8;
    bool hasInfrared = m_b.pointFormatId == 8;

    uint8_t returnNumber(1);
    uint8_t numberOfReturns(1);
    if (loc.returnNumber.present())
        returnNumber = pt.get<uint8_t>(loc.returnNumber);
    if (loc.numberOfReturns.present())
        numberOfReturns = pt.get<uint8_t>(loc.numberOfReturns);

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

    double x = pt.get<double>(loc.x);
    int32_t xi = converter((x - m_b.xform.offset.x) / m_b.xform.scale.x, Dimension::Id::X);
    double y = pt.get<double>(loc.y);
    int32_t yi = converter((y - m_b.xform.offset.y) / m_b.xform.scale.y, Dimension::Id::Y);
    double z = pt.get<double>(loc.z);
    int32_t zi = converter((z - m_b.xform.offset.z) / m_b.xform.scale.z, Dimension::Id::Z);

    ostream << xi << yi << zi;

    ostream << pt.get<uint16_t>(loc.intensity);
    ostream << (uint8_t)(returnNumber | (numberOfReturns << 4));
    ostream << pt.get<uint8_t>(loc.bits);
    ostream << pt.get<uint8_t>(loc.classification);

    uint8_t userData = pt.get<uint8_t>(loc.userData);
    int16_t scanAngleRank =
        static_cast<int16_t>(std::round(pt.get<float>(loc.scanAngleRank) / .006f));
    ostream << userData << scanAngleRank;

    ostream << pt.get<uint16_t>(loc.pointSourceId);

    if (hasTime)
        ostream << pt.get<double>(loc.gpsTime);

    if (hasColor)
    {
        ostream << pt.get<uint16_t>(loc.red);
        ostream << pt.get<uint16_t>(loc.green);
        ostream << pt.get<uint16_t>(loc.blue);
    }

    if (hasInfrared)
        ostream << pt.get<uint16_t>(loc.infrared);

    for (const Field& eloc : extraLocs)
    {
        Everything e;
        std::memset(&e, 0, sizeof(e));
        pt.copyRaw(eloc, &e);
        Utils::insertDim(ostream, eloc.type, e);
    }
}

} // namespace bu
} // namespace untwine
