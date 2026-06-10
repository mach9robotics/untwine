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

#pragma once

#include <memory>

#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>

#include "Stats.hpp"
#include "../untwine/Common.hpp"
#include "../untwine/ThreadPool.hpp"
#include "../untwine/VoxelKey.hpp"

namespace untwine
{
namespace bu
{

class PyramidManager;

constexpr int NumChunkWriters = 8;
constexpr int ChunkWriterQueueSize = 16;

// Compresses octant point data and writes it into the single COPC output file.
//
// A Processor's completion signal unblocks its parent octant, but the compressed
// output it produces is read by nothing until the file is finalized. So rather than
// compress/write on the Processor thread (delaying the parent), Processors enqueue
// their fully-built point views here and signal completion immediately. Worker
// threads compress each chunk, reserve a file location (PyramidManager::newChunk)
// and write it. Chunks land in the file in whatever order they complete - each
// gets a disjoint byte range, so ordering doesn't matter.
//
// The queue is bounded (ChunkWriterQueueSize): each pending chunk pins its point
// data in memory, so enqueue() blocks when workers fall behind. stop() drains the
// queue completely - the chunk table, hierarchy and stats are only complete once
// every chunk has been processed, so it must be called before the file is finalized.
class ChunkWriter
{
public:
    struct Chunk
    {
        VoxelKey key;
        // The view refers to the table, so the table must stay alive until the
        // chunk is written.
        std::shared_ptr<pdal::PointTable> table;
        pdal::PointViewPtr view;
        pdal::DimTypeList extraDims;
        IndexedStats stats;
        size_t count;
    };

    ChunkWriter(PyramidManager& manager, const BaseInfo& b);

    // Blocks if the queue is full.
    void enqueue(Chunk&& chunk);
    // Process all queued chunks and stop the worker threads. Throws if any
    // worker had an error.
    void stop();

private:
    void write(Chunk& chunk);
    void createChunk(const Chunk& chunk);
    void sortChunk(pdal::PointViewPtr view);
    void fillPointBuf(pdal::PointRef& point, std::vector<char>& buf,
        pdal::Dimension::Id bitsDim, const pdal::DimTypeList& extraDims);

    PyramidManager& m_manager;
    const BaseInfo& m_b;
    ThreadPool m_pool;
};

} // namespace bu
} // namespace untwine
