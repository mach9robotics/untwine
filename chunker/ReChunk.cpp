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
#include <memory>
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

// Levels one split step descends: the smallest d with 8^d >= count/target, so a uniform spread
// lands each child at or under the target in a single pass (the trick localGridDepth uses for the
// in-RAM build). Capped by MaxSplitStepDepth to bound the open files and staging buffers, and by
// the headroom left below MaxRechunkLevel.
int splitStepDepth(uint64_t count, uint64_t target, int level)
{
    const uint64_t ratio = (count + target - 1) / target;
    int d = 1;
    while (d < MaxSplitStepDepth && (uint64_t(1) << (3 * d)) < ratio)
        ++d;
    if (d > MaxRechunkLevel - level)
        d = MaxRechunkLevel - level;
    return d < 1 ? 1 : d;
}

// Split `<key>.bin` several levels in ONE read+write pass: each record's descendant cell `d`
// levels down is found by walking the octant midpoints, and the record is dealt directly into
// that cell's .bin. One pass regardless of depth is the point — a per-level split rewrites the
// same bytes at every level, and the first level of a monster chunk is inherently serial, so
// levels are the thing to economize. Returns the children created, each with its point count;
// empty when the file is terminal — within `target`, or at the level cap (a coincident cluster
// never separates; the Indexing phase guards the residual case). Each call touches only its own
// subtree's files, so distinct keys split concurrently.
std::vector<std::pair<VoxelKey, uint64_t>> splitChunkStep(const std::string& tempDir,
    const pdal::BOX3D& fullBounds, size_t pointSize, const VoxelKey& key, uint64_t target)
{
    const std::string path = binPath(tempDir, key);
    const uint64_t count = pdal::FileUtils::fileSize(path) / pointSize;

    if (count <= target || key.level() >= MaxRechunkLevel)
        return {};

    const pdal::BOX3D nodeBounds = voxelBounds(fullBounds, key);
    const int depth = splitStepDepth(count, target, key.level());

    // Cells are created lazily on first point: a clustered chunk occupies few of the 8^d slots.
    // Records stage in per-cell buffers flushed in blocks — a per-record ofstream::write costs
    // more in stream bookkeeping than in I/O.
    struct Cell
    {
        VoxelKey key;
        std::ofstream out;
        std::vector<char> buf;
        uint64_t count = 0;
    };
    std::vector<std::unique_ptr<Cell>> cells(size_t(1) << (3 * depth));

    forEachBlock(path, pointSize, [&](const char *block, uint64_t n)
    {
        for (uint64_t i = 0; i < n; ++i)
        {
            const char *rec = block + i * pointSize;
            const Point p(const_cast<char *>(rec));

            // Walk down `depth` levels of octant midpoints, same subdivision as voxelBounds.
            VoxelKey k = key;
            pdal::BOX3D b = nodeBounds;
            size_t slot = 0;
            for (int lvl = 0; lvl < depth; ++lvl)
            {
                const int oct = bu::childOctant(p.x(), p.y(), p.z(), b);
                const double mx = (b.minx + b.maxx) * 0.5;
                const double my = (b.miny + b.maxy) * 0.5;
                const double mz = (b.minz + b.maxz) * 0.5;
                if (oct & 1) b.minx = mx; else b.maxx = mx;
                if (oct & 2) b.miny = my; else b.maxy = my;
                if (oct & 4) b.minz = mz; else b.maxz = mz;
                k = k.child(oct);
                slot = slot * 8 + oct;
            }

            std::unique_ptr<Cell>& c = cells[slot];
            if (!c)
            {
                c.reset(new Cell);
                c->key = k;
                const std::string cellPath = binPath(tempDir, k);
                c->out.open(os::toNative(cellPath), std::ios::binary | std::ios::trunc);
                if (!c->out)
                    throw FatalError("Couldn't open '" + cellPath + "' for output.");
                c->buf.reserve(SplitStageBytes);
            }
            c->buf.insert(c->buf.end(), rec, rec + pointSize);
            c->count++;
            if (c->buf.size() + pointSize > SplitStageBytes)
            {
                c->out.write(c->buf.data(), c->buf.size());
                c->buf.clear();
            }
        }
    });

    std::vector<std::pair<VoxelKey, uint64_t>> children;
    for (std::unique_ptr<Cell>& c : cells)
    {
        if (!c)
            continue;
        if (!c->buf.empty())
            c->out.write(c->buf.data(), c->buf.size());
        if (!c->out)
            throw FatalError("Write failed re-chunking '" + path + "'.");
        children.emplace_back(c->key, c->count);
    }
    pdal::FileUtils::deleteFile(path);
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
            splitChunkStep(tempDir, fullBounds, pointSize, k, target);
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
                splitChunkStep(tempDir, fullBounds, pointSize, key, target);
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
