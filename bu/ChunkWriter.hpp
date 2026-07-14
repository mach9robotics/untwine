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

#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include <pdal/Dimension.hpp>
#include <pdal/DimType.hpp>

#include "Stats.hpp"
#include "../untwine/Common.hpp"
#include "../untwine/RecordLayout.hpp"
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
    // One built node's point data queued for compression + writing. A self-contained unit of work
    // for a writer thread: it owns its bytes, so it outlives the producer and its mmap'd files.
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
    // The standard-LAS field locations fillPointBuf emits, in point-format order. Resolved once in
    // the ctor from m_layout (absent Field if the dimension isn't stored), so the per-point loop
    // reads through a pre-resolved Field instead of hashing the layout per field per point.
    struct StdLocs
    {
        Field x, y, z, intensity, returnNumber, numberOfReturns, bits, classification,
            userData, scanAngleRank, pointSourceId, gpsTime, red, green, blue, infrared;
    };

    // Worker-thread entry for one chunk: accumulate its per-dimension stats from the packed
    // records, build and write the compressed block (createChunk), then log the node to the manager.
    void write(Chunk& chunk);
    // Build the LAZ-compressed point block for `chunk` (extra-byte sizing, GPS-time sort), reserve
    // its file location via PyramidManager::newChunk, and write it.
    void createChunk(const Chunk& chunk);
    // Resolve m_stdLocs from m_layout. Run-invariant, so called once in the ctor.
    void resolveStdLocs();
    // Emit one packed point `pt` into the LAZ point buffer `buf`: standard fields from `loc`,
    // extra dimensions from `extraLocs`.
    void fillPointBuf(const PackedPoint& pt, std::vector<char>& buf, const StdLocs& loc,
        const std::vector<Field>& extraLocs);

    PyramidManager& m_manager;
    const BaseInfo& m_b;
    ThreadPool m_pool;

    // The packed-record byte layout, built once from m_b.dimInfo. Owns pointSize and the
    // id/name -> Field table; produces the PackedPoint views every field read goes through.
    RecordLayout m_layout;
    // The standard fillPointBuf field locations, resolved once (run-invariant) in the ctor.
    StdLocs m_stdLocs;
};

} // namespace bu
} // namespace untwine
