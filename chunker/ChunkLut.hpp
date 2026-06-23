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

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "../untwine/VoxelKey.hpp"

namespace untwine
{
namespace chunker
{

// The chunk partition computed from the count pass's histogram. `roots` is the set of chunk keys;
// each gets one raw-point .bin from the distribute pass and one ChunkBuilder. `cellToRoot` maps
// every occupied count-grid cell to the chunk its points belong to. Every occupied cell maps to
// exactly one chunk: a coarser merged ancestor where the neighborhood is sparse, or the cell
// itself where it can't merge up, so dense regions chunk at the count-grid level.
struct ChunkLut
{
    std::vector<VoxelKey> roots;
    std::unordered_map<VoxelKey, VoxelKey> cellToRoot;
};

// Build the chunk partition by merging the count-grid histogram under a budget. Reuses
// bu::planChunkRoots for the pyramid merge.
//
// cellCounts: occupied count-grid cell -> its point count, from the count pass.
// target:     per-chunk point budget.
// Returns:    the partition: the chunk roots and the cell->chunk map.
ChunkLut buildChunkLut(const std::unordered_map<VoxelKey, uint64_t>& cellCounts, uint64_t target);

} // namespace chunker
} // namespace untwine
