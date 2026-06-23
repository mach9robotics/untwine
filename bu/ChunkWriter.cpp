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
           envOverride("UNTWINE_CHUNK_QUEUE_SIZE", ChunkWriterQueueSize))
{
    m_pool.trap(true, "Unknown error in ChunkWriter");

    // Build the id -> {offset, type} map for a packed temp record straight from m_b.dimInfo, whose
    // canonical ids are set in prep/FilePrep::fillMetadata. The stats loop reads each stat dim by
    // offset through this, and resolveStdLocs() builds m_stdLocs from it.
    for (const FileDimInfo& fdi : m_b.dimInfo)
        m_locById[fdi.dim] = FieldLoc{fdi.offset, fdi.type};

    resolveStdLocs();
}

// Resolve the standard fillPointBuf field locations from m_b.dimInfo. Run-invariant (depends only
// on the layout, fixed for the whole job), so done once here rather than per chunk. Most come from
// m_locById by id; the packed-bits field has a dynamically-assigned id, so it's matched by name.
void ChunkWriter::resolveStdLocs()
{
    using namespace pdal;

    auto byId = [this](Dimension::Id id) -> FieldLoc
    {
        auto it = m_locById.find(id);
        return it != m_locById.end() ? it->second : FieldLoc{};
    };
    m_stdLocs.x = byId(Dimension::Id::X);
    m_stdLocs.y = byId(Dimension::Id::Y);
    m_stdLocs.z = byId(Dimension::Id::Z);
    m_stdLocs.intensity = byId(Dimension::Id::Intensity);
    m_stdLocs.returnNumber = byId(Dimension::Id::ReturnNumber);
    m_stdLocs.numberOfReturns = byId(Dimension::Id::NumberOfReturns);
    m_stdLocs.classification = byId(Dimension::Id::Classification);
    m_stdLocs.userData = byId(Dimension::Id::UserData);
    m_stdLocs.scanAngleRank = byId(Dimension::Id::ScanAngleRank);
    m_stdLocs.pointSourceId = byId(Dimension::Id::PointSourceId);
    m_stdLocs.gpsTime = byId(Dimension::Id::GpsTime);
    m_stdLocs.red = byId(Dimension::Id::Red);
    m_stdLocs.green = byId(Dimension::Id::Green);
    m_stdLocs.blue = byId(Dimension::Id::Blue);
    m_stdLocs.infrared = byId(Dimension::Id::Infrared);
    for (const FileDimInfo& fdi : m_b.dimInfo)
        if (fdi.name == UntwineBitsDimName)
            m_stdLocs.bits = FieldLoc{fdi.offset, fdi.type};
}

namespace
{

// Read sizeof(Src) bytes at p as Src, then numerically cast to T. memcpy (not a pointer cast)
// handles unaligned offsets and avoids a strict-aliasing violation.
template<typename Src, typename T>
inline T readAs(const char* p)
{
    Src v;
    std::memcpy(&v, p, sizeof(v));
    return static_cast<T>(v);
}

} // unnamed namespace

template<typename T>
T ChunkWriter::readField(const char* rec, int offset, pdal::Dimension::Type stored)
{
    using Type = pdal::Dimension::Type;

    // Absent dimension: mirror PDAL's getFieldAs on a missing dim, which yields a zero value.
    if (offset < 0)
        return T{};

    // The stored type is only known at runtime, and converting its bytes to T is a numeric cast
    // (not a reinterpret), so we must dispatch on it. This is getFieldAs<T> minus the layout lookup.
    const char* p = rec + offset;
    switch (stored)
    {
    case Type::Unsigned8:  return readAs<uint8_t,  T>(p);
    case Type::Unsigned16: return readAs<uint16_t, T>(p);
    case Type::Unsigned32: return readAs<uint32_t, T>(p);
    case Type::Unsigned64: return readAs<uint64_t, T>(p);
    case Type::Signed8:    return readAs<int8_t,   T>(p);
    case Type::Signed16:   return readAs<int16_t,  T>(p);
    case Type::Signed32:   return readAs<int32_t,  T>(p);
    case Type::Signed64:   return readAs<int64_t,  T>(p);
    case Type::Float:      return readAs<float,    T>(p);
    case Type::Double:     return readAs<double,   T>(p);
    case Type::None:       return T{};
    }
    return T{};
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
    // (id, Stats) pair keys on the same id the producer registered, so we look up that id's byte
    // offset/type in m_locById and read the value straight from the packed record.
    if (m_b.opts.stats)
    {
        const size_t pointSize = m_b.pointSize;
        for (auto& sp : chunk.stats)
        {
            pdal::Dimension::Id dim = sp.first;
            Stats& s = sp.second;
            auto it = m_locById.find(dim);
            if (it == m_locById.end() || it->second.offset < 0)
                continue;
            const FieldLoc& loc = it->second;
            for (size_t i = 0; i < chunk.count; ++i)
            {
                const char* rec = chunk.records.data() + i * pointSize;
                s.insert(readField<double>(rec, loc.offset, loc.type));
            }
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

    const size_t pointSize = m_b.pointSize;

    // Extra-byte size from the dim types directly; no PointLayout needed.
    int ebCount {0};
    std::vector<FieldLoc> extraLocs;
    extraLocs.reserve(c.extraDims.size());
    for (const DimType& dim : c.extraDims)
    {
        ebCount += (int)Dimension::size(dim.m_type);
        auto it = m_locById.find(dim.m_id);
        // The packed record stores the extra dim with its own offset/type. Use the stored type
        // from the record for the read; the LAZ extra-byte layout uses dim.m_type (they match).
        int offset = (it != m_locById.end()) ? it->second.offset : -1;
        extraLocs.push_back(FieldLoc{offset, dim.m_type});
    }

    // Sort on GPS time, replacing the old PDAL SortFilter. We compute a stable ascending order of
    // point indices by the GpsTime double read at its packed offset, then emit in that order. If
    // there's no GpsTime dimension, emit in natural order.
    std::vector<uint32_t> order(c.count);
    std::iota(order.begin(), order.end(), 0u);
    if (m_stdLocs.gpsTime.offset >= 0)
    {
        const char* base = c.records.data();
        const int gpsOff = m_stdLocs.gpsTime.offset;
        std::stable_sort(order.begin(), order.end(),
            [base, pointSize, gpsOff](uint32_t a, uint32_t b)
            {
                double ta, tb;
                std::memcpy(&ta, base + a * pointSize + gpsOff, sizeof(double));
                std::memcpy(&tb, base + b * pointSize + gpsOff, sizeof(double));
                return ta < tb;
            });
    }

    std::vector<char> buf(lazperf::baseCount(m_b.pointFormatId) + ebCount);
    lazperf::writer::chunk_compressor compressor(m_b.pointFormatId, ebCount);
    for (uint32_t idx : order)
    {
        const char* rec = c.records.data() + (size_t)idx * pointSize;
        fillPointBuf(rec, buf, m_stdLocs, extraLocs);
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

void ChunkWriter::fillPointBuf(const char* rec, std::vector<char>& buf, const StdLocs& loc,
    const std::vector<FieldLoc>& extraLocs)
{
    using namespace pdal;

    LeInserter ostream(buf.data(), buf.size());

    bool hasTime = true;
    bool hasColor = m_b.pointFormatId == 7 || m_b.pointFormatId == 8;
    bool hasInfrared = m_b.pointFormatId == 8;

    uint8_t returnNumber(1);
    uint8_t numberOfReturns(1);
    if (loc.returnNumber.offset >= 0)
        returnNumber = readField<uint8_t>(rec, loc.returnNumber.offset, loc.returnNumber.type);
    if (loc.numberOfReturns.offset >= 0)
        numberOfReturns =
            readField<uint8_t>(rec, loc.numberOfReturns.offset, loc.numberOfReturns.type);

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

    double x = readField<double>(rec, loc.x.offset, loc.x.type);
    int32_t xi = converter((x - m_b.xform.offset.x) / m_b.xform.scale.x, Dimension::Id::X);
    double y = readField<double>(rec, loc.y.offset, loc.y.type);
    int32_t yi = converter((y - m_b.xform.offset.y) / m_b.xform.scale.y, Dimension::Id::Y);
    double z = readField<double>(rec, loc.z.offset, loc.z.type);
    int32_t zi = converter((z - m_b.xform.offset.z) / m_b.xform.scale.z, Dimension::Id::Z);

    ostream << xi << yi << zi;

    ostream << readField<uint16_t>(rec, loc.intensity.offset, loc.intensity.type);
    ostream << (uint8_t)(returnNumber | (numberOfReturns << 4));
    ostream << readField<uint8_t>(rec, loc.bits.offset, loc.bits.type);
    ostream << readField<uint8_t>(rec, loc.classification.offset, loc.classification.type);

    uint8_t userData = readField<uint8_t>(rec, loc.userData.offset, loc.userData.type);
    int16_t scanAngleRank =
        static_cast<int16_t>(std::round(
            readField<float>(rec, loc.scanAngleRank.offset, loc.scanAngleRank.type) / .006f));
    ostream << userData << scanAngleRank;

    ostream << readField<uint16_t>(rec, loc.pointSourceId.offset, loc.pointSourceId.type);

    if (hasTime)
        ostream << readField<double>(rec, loc.gpsTime.offset, loc.gpsTime.type);

    if (hasColor)
    {
        ostream << readField<uint16_t>(rec, loc.red.offset, loc.red.type);
        ostream << readField<uint16_t>(rec, loc.green.offset, loc.green.type);
        ostream << readField<uint16_t>(rec, loc.blue.offset, loc.blue.type);
    }

    if (hasInfrared)
        ostream << readField<uint16_t>(rec, loc.infrared.offset, loc.infrared.type);

    for (const FieldLoc& eloc : extraLocs)
    {
        Everything e;
        std::memset(&e, 0, sizeof(e));
        if (eloc.offset >= 0)
            std::memcpy(&e, rec + eloc.offset, Dimension::size(eloc.type));
        Utils::insertDim(ostream, eloc.type, e);
    }
}

} // namespace bu
} // namespace untwine
