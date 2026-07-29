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

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include <pdal/util/Bounds.hpp>

#include "../untwine/VoxelKey.hpp"

namespace untwine
{
namespace chunker
{

struct ChunkLut;

// Chunking step 3, re-chunk: the plan's per-chunk budget is soft — a single count-grid cell over
// the budget becomes its own chunk whole (planChunkRoots' oversized-leaf case and buildChunkLut's
// fallback), so a dense or coincident-heavy region can produce a chunk far larger than the budget.
// The Indexing phase loads a whole chunk into RAM and indexes it with 32-bit point indices, so an
// unbounded chunk is both a thrash risk and, past INT_MAX points, undefined behavior. These
// functions restore the budget after the distribute pass by recursively splitting any oversized
// chunk .bin into its child octants — a packed-point read/route/append with no PDAL decode,
// mirroring the EPF Reprocessor. A truly coincident cluster can never split, so recursion also
// stops at MaxRechunkLevel; the Indexing phase guards the residual case.

// The bounds of voxel `key` within `fullBounds`: the same subdivision bu::VoxelInfo::bounds()
// computes, without dragging in that class's sampling machinery.
pdal::BOX3D voxelBounds(const pdal::BOX3D& fullBounds, const VoxelKey& key);

// One sequential pass over a chunk .bin: how many of its `count` packed points fall in each child
// octant of `nodeBounds`. Octant order matches VoxelKey::child / bu::childOctant.
std::array<uint64_t, 8> countOctants(const char *data, uint64_t count, size_t pointSize,
    const pdal::BOX3D& nodeBounds);

// Recursively split `<key>.bin` in `tempDir` until every resulting chunk holds at most `target`
// points or sits at MaxRechunkLevel. A split with a single occupied octant is a rename, not a
// rewrite, so a tight cluster descends cheaply until it either separates or hits the level cap.
// Returns the keys of the resulting chunk .bins (just `key` if no split was needed).
std::vector<VoxelKey> rechunkFile(const std::string& tempDir, const pdal::BOX3D& fullBounds,
    size_t pointSize, const VoxelKey& key, uint64_t target);

// Split every chunk whose planned point count exceeds `trigger` down to at most `target` points
// each. The thresholds are deliberately distinct: `trigger` (RechunkTriggerPoints) says when a
// chunk is dangerous, `target` (the planning budget) says what to split it into once the rewrite
// is being paid for anyway. Chunk counts are derived from the count-pass histogram and the LUT,
// so untouched chunks cost nothing. Returns the number of chunks that were split.
size_t rechunkOversized(const std::string& tempDir, const pdal::BOX3D& fullBounds,
    size_t pointSize, const ChunkLut& lut,
    const std::unordered_map<VoxelKey, uint64_t>& cellCounts, uint64_t trigger, uint64_t target);

} // namespace chunker
} // namespace untwine
