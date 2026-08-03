/*****************************************************************************
 *   Copyright (c) 2020, Hobu, Inc. (info@hobu.co)                           *
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

#include <vector>

namespace untwine
{

struct BaseInfo;
struct FileInfo;
class ProgressWriter;
namespace epf { class Grid; }

namespace chunker
{

struct ChunkLut;

// Distribute pass, the second decode: re-read every input point, route it through the count grid
// and LUT to its chunk, and append its full packed record to that chunk's <chunkRoot>.bin via the
// EPF buffered writer. Produces one raw-point .bin per chunk for bu::indexChunks to consume.
//
// b:         shared pipeline state.
// fileInfos: input files to decode. Reordered largest-first as a side effect.
// grid:      the count grid, identical to the one the count pass used.
// lut:       the cell->chunk map from buildChunkLut.
// progress:  progress writer for the 20-40% distribute band.
void distribute(BaseInfo& b, std::vector<FileInfo>& fileInfos, const epf::Grid& grid,
    const ChunkLut& lut, ProgressWriter& progress);

} // namespace chunker
} // namespace untwine
