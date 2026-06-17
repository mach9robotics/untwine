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

#include "Constants.hpp"
#include "../untwine/VoxelKey.hpp"

namespace untwine
{

struct FileInfo;
struct Transform;
class ProgressWriter;
namespace epf { class Grid; }

namespace chunker
{

// Count-grid resolution scaled by total point count: under 100 M points uses level 7 = 128^3,
// under 500 M level 8 = 256^3, otherwise level 9 = 512^3. A finer grid for larger clouds keeps
// occupied-cell counts, and thus the chunk plan, useful.
//
// numPoints: total point count across all input files.
// Returns:   the count-grid octree level.
inline int countGridLevel(uint64_t numPoints)
{
    if (numPoints < CountGridSmallMaxPoints)
        return CountGridLevelSmall;
    if (numPoints < CountGridMediumMaxPoints)
        return CountGridLevelMedium;
    return CountGridLevelLarge;
}

// Count pass: decode every input point once, quantize it to the output transform, and tally how
// many fall in each cell of `grid`. No .bin files are written; this is purely the histogram the
// chunk planner needs.
//
// fileInfos:  input files to decode. Not modified.
// grid:       the count grid whose cell each point is tallied into.
// xform:      output transform the points are quantized to before keying. Must match the
//             distribute pass exactly, or points land in different cells across the two passes.
// numThreads: number of parallel decode threads.
// progress:   ticked once per ChunkSize points so the count decode stays visible to the profiler.
//             The caller scales the incrementer to give it its share of the front-end band.
// Returns:    a map from each occupied count-grid cell to its point count.
std::unordered_map<VoxelKey, uint64_t> countPoints(std::vector<FileInfo>& fileInfos,
    const epf::Grid& grid, const Transform& xform, std::size_t numThreads,
    ProgressWriter& progress);

} // namespace chunker
} // namespace untwine
