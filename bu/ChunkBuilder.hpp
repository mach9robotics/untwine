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

// Chunk-local in-RAM bottom-up build for the --chunked_build path.
//
// The default BU schedules one Processor per octree node bottom-up across the whole tree,
// with an 8-children barrier. Near the root few nodes are independent, so it is span-bound.
//
// runChunkedBuild() instead partitions the leaves into adaptively-sized, spatially-disjoint
// chunks, each a subtree rooted at a level chosen so it holds <= maxChunkPoints points, and
// processes each chunk independently in parallel:
//   * load the chunk's descendant leaf .bins into one PointAccessor by mmap,
//   * sample the subtree bottom-up in RAM, passing accepted points up as index lists with no
//     interior .bin writes and no global barrier,
//   * emit a COPC chunk via PyramidManager for every node strictly below the chunk root,
//   * write the chunk root's accepted points to <root>.bin,
//   * delete the consumed leaf .bins.
//
// After it returns, the temp dir holds only chunk-root .bins, both internal roots and any
// single-leaf roots, which the caller hands to the existing scheduler to merge into the global
// octree from chunk-root level up to the root. One PyramidManager and CopcSupport is shared across
// both phases, so chunk-table, hierarchy, stats, and finalization are unified.
//
// b:             shared pipeline state.
// mgr:           pyramid manager that records each emitted chunk and drives the merge.
// leaves:        leaf VoxelKey mapped to its .bin, from BuPyramid::getInputFiles().
// maxChunkPoints: per-chunk point budget, where 0 lets the caller pick a default.
// progress:      progress writer for the build band, or null.
void runChunkedBuild(const BaseInfo& b, PyramidManager& mgr,
    const std::unordered_map<VoxelKey, FileInfo>& leaves,
    uint64_t maxChunkPoints, ProgressWriter* progress);

// Indexing phase for the chunker pipeline: each entry of `chunkFiles` is already one chunk, a
// <chunkRoot>.bin of raw points produced by Chunking via chunker::distribute. Build a local octree
// from each in parallel, leaving the chunk-root .bins for the Merging phase. Unlike
// runChunkedBuild, there is no re-planning, since Chunking already decided the partition.
//
// b:          shared pipeline state.
// mgr:        pyramid manager that records each emitted chunk and drives the merge.
// chunkFiles: chunk-root VoxelKey mapped to its raw-point .bin, from Chunking.
// progress:   progress writer for the indexing band, or null.
void indexChunks(const BaseInfo& b, PyramidManager& mgr,
    const std::unordered_map<VoxelKey, FileInfo>& chunkFiles, ProgressWriter* progress);

} // namespace bu
} // namespace untwine
