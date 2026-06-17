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

#include <vector>

namespace untwine
{

struct BaseInfo;
struct FileInfo;
class ProgressWriter;

namespace chunker
{

// Chunking phase of the counting-sort pipeline: count, plan, then distribute. Count tallies points
// into the fine count grid; plan merges those cells bottom-up under a per-chunk budget so sparse
// neighborhoods collapse into one chunk root, yielding far fewer chunks than occupied cells;
// distribute re-reads the points and writes one raw-point .bin per chunk for the Indexing phase.
// An alternative to EPF binning and the iterative Reprocessor. Gated behind --chunker so it A/B's
// against EPF.
class Chunker
{
public:
    // common: shared pipeline state: options, bounds, dims, and output transform.
    Chunker(BaseInfo& common) : m_b(common)
    {}

    // Run the chunking phase: count, plan, then distribute.
    //
    // progress:  progress writer for the 0-40% chunking band.
    // fileInfos: input files to chunk. Reordered largest-first as a side effect.
    void run(ProgressWriter& progress, std::vector<FileInfo>& fileInfos);

private:
    BaseInfo& m_b;
};

} // namespace chunker
} // namespace untwine
