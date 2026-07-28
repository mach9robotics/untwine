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

#include <cmath>
#include <fstream>

#include <pdal/util/FileUtils.hpp>

#include "ReChunk.hpp"
#include "ChunkLut.hpp"
#include "Constants.hpp"

#include "../untwine/Common.hpp"
#include "../untwine/Point.hpp"
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

std::vector<VoxelKey> rechunkFile(const std::string& tempDir, const pdal::BOX3D& fullBounds,
    size_t pointSize, const VoxelKey& key, uint64_t target)
{
    const std::string path = binPath(tempDir, key);
    const uint64_t count = pdal::FileUtils::fileSize(path) / pointSize;

    // Under budget, or too deep to keep splitting (a coincident cluster never separates): this is
    // a final chunk. The Indexing phase splits what it can in RAM and guards the residual case.
    if (count <= target || key.level() >= MaxRechunkLevel)
        return { key };

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
        return rechunkFile(tempDir, fullBounds, pointSize, child, target);
    }

    // Pass 2: route each record to its octant's .bin. Occupancy is known, so the output streams
    // can be opened up front; ofstream's own buffering batches the per-record appends.
    {
        std::array<std::ofstream, 8> outs;
        for (int i = 0; i < 8; ++i)
        {
            if (!counts[i])
                continue;
            const std::string childPath = binPath(tempDir, key.child(i));
            outs[i].open(os::toNative(childPath), std::ios::binary | std::ios::trunc);
            if (!outs[i])
                throw FatalError("Couldn't open '" + childPath + "' for output.");
        }
        forEachBlock(path, pointSize, [&](const char *block, uint64_t n)
        {
            for (uint64_t i = 0; i < n; ++i)
            {
                const char *rec = block + i * pointSize;
                const Point p(const_cast<char *>(rec));
                outs[bu::childOctant(p.x(), p.y(), p.z(), nodeBounds)].write(rec, pointSize);
            }
        });
        for (int i = 0; i < 8; ++i)
            if (counts[i] && !outs[i])
                throw FatalError("Write failed re-chunking '" + path + "'.");
    }
    pdal::FileUtils::deleteFile(path);

    std::vector<VoxelKey> chunks;
    for (int i = 0; i < 8; ++i)
    {
        if (!counts[i])
            continue;
        std::vector<VoxelKey> sub =
            rechunkFile(tempDir, fullBounds, pointSize, key.child(i), target);
        chunks.insert(chunks.end(), sub.begin(), sub.end());
    }
    return chunks;
}

size_t rechunkOversized(const std::string& tempDir, const pdal::BOX3D& fullBounds,
    size_t pointSize, const ChunkLut& lut,
    const std::unordered_map<VoxelKey, uint64_t>& cellCounts, uint64_t target)
{
    // Planned point count per chunk. Merged (internal) roots hold <= target by construction, so
    // only oversized count-grid cells that became their own chunks can exceed it.
    std::unordered_map<VoxelKey, uint64_t> rootCounts;
    for (const auto& c : cellCounts)
        rootCounts[lut.cellToRoot.at(c.first)] += c.second;

    size_t split = 0;
    for (const auto& rc : rootCounts)
    {
        if (rc.second <= target)
            continue;
        rechunkFile(tempDir, fullBounds, pointSize, rc.first, target);
        ++split;
    }
    return split;
}

} // namespace chunker
} // namespace untwine
