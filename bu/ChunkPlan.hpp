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

#include <cstdint>
#include <unordered_map>
#include <vector>

#include <pdal/util/Bounds.hpp>

#include "../untwine/VoxelKey.hpp"

namespace untwine
{
namespace bu
{

// Pick the internal chunk roots for the chunk-local build. Descends from the global root and
// stops at the shallowest key whose subtree holds <= `target` points. That key becomes a chunk
// root unless it's a leaf; a single-leaf chunk root is handled directly by the merge and is not
// returned. The global root is never returned, so chunk roots are at level >= 1. Pure logic with
// no I/O or threading deps. Unit tested in test/ChunkBuildTest.cpp.
//
// leafCounts: each leaf VoxelKey mapped to its point count.
// target:     per-chunk point budget.
// Returns:    the chunk-root keys.
std::vector<VoxelKey> planChunkRoots(const std::unordered_map<VoxelKey, uint64_t>& leafCounts,
    uint64_t target);

// Depth in levels below a chunk root for the per-chunk local count grid that the chunk-local build
// histograms points into before merging cells under the leaf budget. Chosen so a uniform spread
// puts about `leafBudget` points in each finest cell, plus one extra level of margin. Cells denser
// than the budget are refined by recursing with a fresh grid, so this is the per-step granularity,
// not a hard resolution limit. Pure logic, unit tested.
//
// numPoints:         number of points in the chunk being structured.
// leafBudget:        target maximum points per finest cell.
// maxDepthBelowRoot: cap on the returned depth, keeping root depth within the octree limit.
// Returns:           the grid depth, in levels below the chunk root.
int localGridDepth(uint64_t numPoints, uint64_t leafBudget, int maxDepthBelowRoot);

// The child octant of node `b` containing a point, split at the node midpoint. The bit layout
// matches VoxelKey::child: bit 0 is x, bit 1 is y, bit 2 is z. A point exactly on a midpoint goes
// to the upper side, matching how a leaf-boundary point is binned.
//
// x: the point's x coordinate.
// y: the point's y coordinate.
// z: the point's z coordinate.
// b: the node bounds, split at their midpoint.
// Returns: the child octant, 0..7, containing the point.
inline int childOctant(double x, double y, double z, const pdal::BOX3D& b)
{
    const double mx = (b.minx + b.maxx) * 0.5;
    const double my = (b.miny + b.maxy) * 0.5;
    const double mz = (b.minz + b.maxz) * 0.5;
    return (x >= mx ? 1 : 0) | (y >= my ? 2 : 0) | (z >= mz ? 4 : 0);
}

} // namespace bu
} // namespace untwine
