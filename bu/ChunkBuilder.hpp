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

// Indexing phase for the chunker pipeline: each entry of `chunkFiles` is already one chunk, a
// <chunkRoot>.bin of raw points produced by Chunking via chunker::distribute. Build a local octree
// from each in parallel, leaving the chunk-root .bins for the Merging phase. Chunking already
// decided the partition, so there is no re-planning here.
//
// b:          shared pipeline state.
// mgr:        pyramid manager that records each emitted chunk and drives the merge.
// chunkFiles: chunk-root VoxelKey mapped to its raw-point .bin, from Chunking.
// progress:   progress writer for the indexing band, or null.
void indexChunks(const BaseInfo& b, PyramidManager& mgr,
    const std::unordered_map<VoxelKey, FileInfo>& chunkFiles, ProgressWriter* progress);

} // namespace bu
} // namespace untwine
