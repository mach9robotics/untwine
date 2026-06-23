/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                      *
 *                                                                           *
 *   All rights reserved.                                                    *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 ****************************************************************************/

#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include <pdal/Dimension.hpp>
#include <pdal/DimType.hpp>

#include "Stats.hpp"
#include "../untwine/Common.hpp"
#include "../untwine/ThreadPool.hpp"
#include "../untwine/VoxelKey.hpp"

namespace untwine
{
namespace bu
{

class PyramidManager;

// Defaults; overridable at runtime via UNTWINE_NUM_CHUNK_WRITERS / UNTWINE_CHUNK_QUEUE_SIZE
// (see ChunkWriter's constructor and untwine::envOverride).
//
// BU is compression-throughput-bound through the bounded queue, and the writer threads are
// I/O/serialization-bound (per-chunk file write + newChunk mutex), not CPU-bound — so the pool
// is sized to the core count, not below it. Tune via env on other-sized boxes.
constexpr int NumChunkWriters = 16;
constexpr int ChunkWriterQueueSize = 16;

// Compresses octant point data and writes it into the single COPC output file.
//
// Terminology: the "chunks" here are COPC chunks. One octree node's LAZ-compressed block, the unit
// the COPC chunk table indexes. This is unrelated to a chunker/ Chunking-phase chunk (a spatial
// subtree written as one .bin); see chunker/Chunker. @Kyle: we should probably just rename this
// to COPCWriter in the future.
//
// A node's completion signal unblocks its parent, but the compressed output it produces is read
// by nothing until the file is finalized. So rather than compress and write on the build thread
// (delaying the parent), builders enqueue their fully-built point data here and signal completion
// immediately. Worker threads compress each chunk, reserve a file location via
// PyramidManager::newChunk, and write it. Chunks land in the file in whatever order they complete;
// each gets a disjoint byte range, so ordering doesn't matter.
//
// The queue is bounded (ChunkWriterQueueSize): each pending chunk pins its point data in memory,
// so enqueue() blocks when workers fall behind. stop() drains the queue completely; the chunk
// table, hierarchy, and stats are only complete once every chunk has been processed, so it must be
// called before the file is finalized.
//
// Point data is carried as the producer's own packed temp-record bytes (m_b.pointSize per point),
// not a pdal::PointView. The producer bulk-memcpys each point's record into `records`; this writer
// reads fields back out by their known byte offset (built once in the ctor from m_b.dimInfo). This
// avoids ~15 per-field setField/getFieldAs calls per point on both sides; the byte buffer is fully
// self-contained, so deferred compression after the producer's mmap is gone is still safe.
class ChunkWriter
{
public:
    struct Chunk
    {
        VoxelKey key;
        // Packed point records: count * m_b.pointSize bytes, point i at records[i*pointSize].
        // Self-contained (owns its bytes), so it outlives the producer and its mmap'd files.
        std::vector<char> records;
        pdal::DimTypeList extraDims;
        IndexedStats stats;
        size_t count;
    };

    ChunkWriter(PyramidManager& manager, const BaseInfo& b);

    // Enqueue a chunk for compression and writing. Blocks if the queue is full.
    //
    // chunk: the fully-built, self-contained chunk to write (owns its packed records).
    void enqueue(Chunk&& chunk);

    // Process all queued chunks and stop the worker threads. Throws if any worker had an error.
    void stop();

private:
    // The byte offset and stored type of a dimension within a packed temp record.
    struct FieldLoc
    {
        int offset {-1};                  // -1 == dimension absent from the packed record
        pdal::Dimension::Type type {pdal::Dimension::Type::None};
    };

    // The standard-LAS field locations fillPointBuf emits, in point-format order. Resolved once
    // per chunk in createChunk (offset -1 if absent), so the per-point loop reads a plain offset
    // instead of hashing m_locById per field per point.
    struct StdLocs
    {
        FieldLoc x, y, z, intensity, returnNumber, numberOfReturns, bits, classification,
            userData, scanAngleRank, pointSourceId, gpsTime, red, green, blue, infrared;
    };

    void write(Chunk& chunk);
    void createChunk(const Chunk& chunk);
    // Resolve m_stdLocs from m_b.dimInfo. Run-invariant, so called once in the ctor.
    void resolveStdLocs();
    void fillPointBuf(const char* rec, std::vector<char>& buf, const StdLocs& loc,
        const std::vector<FieldLoc>& extraLocs);

    // Read the value stored at rec+offset (interpreted as `stored`) and cast it to T. Robust to
    // a stored type differing from the requested type.
    template<typename T>
    static T readField(const char* rec, int offset, pdal::Dimension::Type stored);

    PyramidManager& m_manager;
    const BaseInfo& m_b;
    ThreadPool m_pool;

    // id -> {offset, type} within a packed record, built once from m_b.dimInfo. Used by the stats
    // loop (keyed on the same pdal::Dimension::Id the producer's IndexedStats use) and to build
    // m_stdLocs.
    std::unordered_map<pdal::Dimension::Id, FieldLoc> m_locById;
    // The standard fillPointBuf field locations, resolved once (run-invariant) in the ctor.
    StdLocs m_stdLocs;
};

} // namespace bu
} // namespace untwine
