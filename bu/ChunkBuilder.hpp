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

#include "../untwine/VoxelKey.hpp"
#include "FileInfo.hpp"

namespace untwine
{

struct BaseInfo;
class ProgressWriter;

namespace bu
{

class PyramidManager;

// Chunk-local in-RAM bottom-up build (the "--chunked-build" path).
//
// The default BU schedules one Processor per octree node bottom-up across the whole tree,
// with an 8-children barrier. Near the root few nodes are independent, so it is span-bound.
//
// runChunkedBuild() instead partitions the leaves into adaptively-sized, spatially-disjoint
// "chunks" (subtrees rooted at a level chosen so each subtree holds <= maxChunkPoints points)
// and processes each chunk independently in parallel:
//   * load the chunk's descendant leaf .bins into one PointAccessor (mmap),
//   * sample the subtree bottom-up *in RAM*, passing accepted points up as index lists
//     (no interior .bin writes, no global barrier),
//   * emit a COPC chunk for every node strictly below the chunk root (via PyramidManager),
//   * write the chunk root's accepted points to <root>.bin,
//   * delete the consumed leaf .bins.
//
// After it returns, the temp dir holds only chunk-root .bins (internal roots plus any
// single-leaf roots), which the caller hands to the existing scheduler to build the small
// "cap" (chunk-root level up to the global root). One PyramidManager/CopcSupport is shared
// across both phases, so chunk-table / hierarchy / stats / finalization are unified.
//
// `leaves` is the leaf map discovered by BuPyramid::getInputFiles(). maxChunkPoints is the
// per-chunk point budget (0 lets the caller pick a default).
void runChunkedBuild(const BaseInfo& b, PyramidManager& mgr,
    const std::unordered_map<VoxelKey, FileInfo>& leaves,
    uint64_t maxChunkPoints, ProgressWriter* progress);

} // namespace bu
} // namespace untwine
