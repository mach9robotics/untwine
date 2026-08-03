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

#include <unordered_set>

#include "ChunkLut.hpp"
#include "../bu/ChunkPlan.hpp"   // bu::planChunkRoots

namespace untwine
{
namespace chunker
{

ChunkLut buildChunkLut(const std::unordered_map<VoxelKey, uint64_t>& cellCounts, uint64_t target)
{
    const VoxelKey root;

    ChunkLut lut;
    lut.roots = bu::planChunkRoots(cellCounts, target);
    std::unordered_set<VoxelKey> rootSet(lut.roots.begin(), lut.roots.end());

    lut.cellToRoot.reserve(cellCounts.size());
    for (const auto& c : cellCounts)
    {
        // Walk up to the nearest chunk root that covers this cell.
        VoxelKey k = c.first;
        VoxelKey chunk;
        bool found = false;
        while (true)
        {
            if (rootSet.count(k))
            {
                chunk = k;
                found = true;
                break;
            }
            if (k == root)
                break;
            k = k.parent();
        }
        // No covering root, so a small cell whose over-budget neighborhood didn't merge it up
        // becomes its own chunk at the count-grid level.
        if (!found)
        {
            chunk = c.first;
            if (rootSet.insert(chunk).second)
                lut.roots.push_back(chunk);
        }
        lut.cellToRoot[c.first] = chunk;
    }
    return lut;
}

} // namespace chunker
} // namespace untwine
