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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <thread>

#include <pdal/util/FileUtils.hpp>

#include "ReChunk.hpp"
#include "ChunkLut.hpp"
#include "Constants.hpp"

#include "../untwine/Common.hpp"
#include "../untwine/Point.hpp"
#include "../untwine/ThreadPool.hpp"
#include "../bu/ChunkPlan.hpp"  // bu::childOctant

#include <stringconv.hpp>  // untwine/os

namespace untwine
{
namespace chunker
{

namespace
{

// Read granularity for the streaming passes. Reads are rounded down to whole point records so a
// record never straddles two reads.
constexpr size_t ReadBlockBytes = 4 * 1024 * 1024;

std::string binPath(const std::string& tempDir, const VoxelKey& key)
{
    return tempDir + "/" + key.toString() + ".bin";
}

// Stream a chunk .bin through `consume(block, count)` one whole-record block at a time.
template <typename Consumer>
void forEachBlock(const std::string& path, size_t pointSize, Consumer consume)
{
    std::ifstream in(os::toNative(path), std::ios::binary);
    if (!in)
        throw FatalError("Couldn't open '" + path + "' for re-chunking.");

    const size_t blockPoints = (std::max)(ReadBlockBytes / pointSize, size_t(1));
    std::vector<char> block(blockPoints * pointSize);
    while (in)
    {
        in.read(block.data(), block.size());
        const size_t bytes = (size_t)in.gcount();
        if (bytes % pointSize)
            throw FatalError("Chunk file '" + path + "' size isn't a multiple of the "
                "point size.");
        if (bytes)
            consume(block.data(), bytes / pointSize);
    }
}

} // unnamed namespace

pdal::BOX3D voxelBounds(const pdal::BOX3D& fullBounds, const VoxelKey& key)
{
    // Same subdivision as bu::VoxelInfo::bounds(): the full bounds divided into 2^level cells per
    // axis, with this voxel's cell selected by the key index.
    const double cells = std::pow(2, key.level());
    const double xw = (fullBounds.maxx - fullBounds.minx) / cells;
    const double yw = (fullBounds.maxy - fullBounds.miny) / cells;
    const double zw = (fullBounds.maxz - fullBounds.minz) / cells;

    pdal::BOX3D b;
    b.minx = fullBounds.minx + key.x() * xw;
    b.maxx = b.minx + xw;
    b.miny = fullBounds.miny + key.y() * yw;
    b.maxy = b.miny + yw;
    b.minz = fullBounds.minz + key.z() * zw;
    b.maxz = b.minz + zw;
    return b;
}

std::array<uint64_t, 8> countOctants(const char *data, uint64_t count, size_t pointSize,
    const pdal::BOX3D& nodeBounds)
{
    std::array<uint64_t, 8> counts{};
    for (uint64_t i = 0; i < count; ++i)
    {
        const Point p(const_cast<char *>(data + i * pointSize));
        counts[bu::childOctant(p.x(), p.y(), p.z(), nodeBounds)]++;
    }
    return counts;
}

// Split `<key>.bin` one level: route its records into the child octants' .bins (or rename when a
// single octant holds everything) and return the children created, each with its point count.
// Returns an empty vector when the file is terminal — already within `target` or at the level
// cap (a coincident cluster never separates; the Indexing phase guards the residual case).
// Each call touches only its own subtree's files, so distinct keys split concurrently.
std::vector<std::pair<VoxelKey, uint64_t>> splitChunkOnce(const std::string& tempDir,
    const pdal::BOX3D& fullBounds, size_t pointSize, const VoxelKey& key, uint64_t target)
{
    const std::string path = binPath(tempDir, key);
    const uint64_t count = pdal::FileUtils::fileSize(path) / pointSize;

    if (count <= target || key.level() >= MaxRechunkLevel)
        return {};

    const pdal::BOX3D nodeBounds = voxelBounds(fullBounds, key);

    // Pass 1: octant occupancy, so a no-op split (everything in one octant) is a rename rather
    // than a rewrite.
    std::array<uint64_t, 8> counts{};
    forEachBlock(path, pointSize, [&](const char *block, uint64_t n)
    {
        const std::array<uint64_t, 8> c = countOctants(block, n, pointSize, nodeBounds);
        for (int i = 0; i < 8; ++i)
            counts[i] += c[i];
    });

    int occupied = 0;
    for (int i = 0; i < 8; ++i)
        if (counts[i])
            ++occupied;
    if (occupied == 1)
    {
        int dir = 0;
        while (!counts[dir])
            ++dir;
        const VoxelKey child = key.child(dir);
        pdal::FileUtils::renameFile(binPath(tempDir, child), path);
        return { { child, count } };
    }

    // Pass 2: route each record to its octant's .bin. Occupancy is known, so the output streams
    // can be opened up front. Records are staged in per-octant block buffers and written out in
    // block-sized chunks — a per-record ofstream::write of a few dozen bytes costs more in stream
    // bookkeeping than in I/O and dominates the split otherwise.
    {
        std::array<std::ofstream, 8> outs;
        std::array<std::vector<char>, 8> bufs;
        for (int i = 0; i < 8; ++i)
        {
            if (!counts[i])
                continue;
            const std::string childPath = binPath(tempDir, key.child(i));
            outs[i].open(os::toNative(childPath), std::ios::binary | std::ios::trunc);
            if (!outs[i])
                throw FatalError("Couldn't open '" + childPath + "' for output.");
            bufs[i].reserve(ReadBlockBytes);
        }
        forEachBlock(path, pointSize, [&](const char *block, uint64_t n)
        {
            for (uint64_t i = 0; i < n; ++i)
            {
                const char *rec = block + i * pointSize;
                const Point p(const_cast<char *>(rec));
                const int oct = bu::childOctant(p.x(), p.y(), p.z(), nodeBounds);
                std::vector<char>& buf = bufs[oct];
                buf.insert(buf.end(), rec, rec + pointSize);
                if (buf.size() + pointSize > ReadBlockBytes)
                {
                    outs[oct].write(buf.data(), buf.size());
                    buf.clear();
                }
            }
        });
        for (int i = 0; i < 8; ++i)
        {
            if (!bufs[i].empty())
                outs[i].write(bufs[i].data(), bufs[i].size());
            if (counts[i] && !outs[i])
                throw FatalError("Write failed re-chunking '" + path + "'.");
        }
    }
    pdal::FileUtils::deleteFile(path);

    std::vector<std::pair<VoxelKey, uint64_t>> children;
    for (int i = 0; i < 8; ++i)
        if (counts[i])
            children.emplace_back(key.child(i), counts[i]);
    return children;
}

std::vector<VoxelKey> rechunkFile(const std::string& tempDir, const pdal::BOX3D& fullBounds,
    size_t pointSize, const VoxelKey& key, uint64_t target)
{
    std::vector<VoxelKey> finals;
    std::vector<VoxelKey> work{ key };
    while (!work.empty())
    {
        const VoxelKey k = work.back();
        work.pop_back();
        std::vector<std::pair<VoxelKey, uint64_t>> children =
            splitChunkOnce(tempDir, fullBounds, pointSize, k, target);
        if (children.empty())
        {
            finals.push_back(k);
            continue;
        }
        for (const auto& c : children)
        {
            if (c.second > target)
                work.push_back(c.first);
            else
                finals.push_back(c.first);
        }
    }
    return finals;
}

size_t rechunkOversized(const std::string& tempDir, const pdal::BOX3D& fullBounds,
    size_t pointSize, const ChunkLut& lut,
    const std::unordered_map<VoxelKey, uint64_t>& cellCounts, uint64_t trigger, uint64_t target)
{
    // A split target above the trigger would recreate over-trigger chunks; clamp.
    target = (std::min)(target, trigger);

    // Planned point count per chunk. Merged (internal) roots hold <= target by construction, so
    // only oversized count-grid cells that became their own chunks can exceed it.
    std::unordered_map<VoxelKey, uint64_t> rootCounts;
    for (const auto& c : cellCounts)
        rootCounts[lut.cellToRoot.at(c.first)] += c.second;

    std::vector<VoxelKey> oversized;
    for (const auto& rc : rootCounts)
        if (rc.second > trigger)
            oversized.push_back(rc.first);
    if (oversized.empty())
        return 0;

    // One task per SPLIT STEP, not per oversized chunk: a task splits its file one level and
    // re-enqueues each still-oversized child as its own task, so the fan-out grows 8x per level
    // and a single monster chunk spreads across the pool instead of pinning one thread (the
    // sparse-bounds case concentrates most of the input in a handful of chunks). Tasks touch
    // only their own subtree's files, so steps are independent. await() is documented to allow
    // add() while awaiting, and returns only when the queue is empty AND nothing is running —
    // a running task enqueues its children before it completes, so no work is missed.
    std::size_t nthreads = std::thread::hardware_concurrency();
    if (nthreads == 0)
        nthreads = DefaultWorkerThreads;
    ThreadPool pool(nthreads);
    pool.trap(true);
    std::function<void(const VoxelKey&)> submit = [&](const VoxelKey& key)
    {
        pool.add([&, key]()
        {
            const std::vector<std::pair<VoxelKey, uint64_t>> children =
                splitChunkOnce(tempDir, fullBounds, pointSize, key, target);
            for (const auto& c : children)
                if (c.second > target)
                    submit(c.first);
        });
    };
    for (const VoxelKey& key : oversized)
        submit(key);
    pool.await();
    if (pool.hasErrors())
    {
        std::vector<std::string> errs = pool.clearErrors();
        throw FatalError(errs.empty() ? std::string("Re-chunk failed.") : errs.front());
    }
    pool.join();
    return oversized.size();
}

} // namespace chunker
} // namespace untwine
