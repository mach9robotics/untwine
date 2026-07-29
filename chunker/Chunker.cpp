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

#include <iostream>
#include <thread>

#include "Chunker.hpp"
#include "ChunkLut.hpp"
#include "Constants.hpp"
#include "CountPass.hpp"
#include "Distribute.hpp"
#include "ReChunk.hpp"

#include "../untwine/Common.hpp"
#include "../untwine/FileInfo.hpp"
#include "../untwine/ProgressWriter.hpp"
#include "pointio/Grid.hpp"

namespace untwine
{
namespace chunker
{

void Chunker::run(ProgressWriter& progress, std::vector<FileInfo>& fileInfos)
{
    const int countLevel = countGridLevel(m_b.numPoints);
    epf::Grid grid(m_b.bounds, m_b.numPoints, countLevel, m_b.opts.doCube);

    std::size_t nthreads = std::thread::hardware_concurrency();
    if (nthreads == 0)
        nthreads = DefaultWorkerThreads;

    // Chunking is two full decodes, count then distribute. Size the 0-40% band over 2 * numPoints
    // so the count step drives 0-20% and distribute 20-40%. Otherwise the count decode emits no
    // progress and is invisible to the profiler's phase windows.
    progress.setPointIncrementer(2 * m_b.numPoints, 40);

    // Chunking step 1, count: decode, quantize, and tally into the count grid. No .bin written.
    std::unordered_map<VoxelKey, uint64_t> cellCounts =
        countPoints(fileInfos, grid, m_b.xform, nthreads, progress);

    uint64_t counted = 0;
    for (const auto& c : cellCounts)
        counted += c.second;

    // Chunking step 2a, plan: pyramid-merge the count histogram into chunk roots and a cell->chunk
    // LUT. This is the same merge-under-count the Indexing phase uses locally, here at chunk scale.
    uint64_t target = m_b.opts.maxChunkPoints ? m_b.opts.maxChunkPoints : DefaultMaxChunkPoints;
    ChunkLut lut = buildChunkLut(cellCounts, target);

    if (m_b.opts.progressDebug)
        std::cerr << "[chunking] count-level " << countLevel
                  << " | occupied cells " << cellCounts.size()
                  << " | counted points " << counted << " (input " << m_b.numPoints << ")"
                  << " | target " << target
                  << " | chunks " << lut.roots.size() << "\n";

    // Chunking step 2b, distribute: decode again and route each point through the LUT into its
    // chunk's .bin. The Indexing phase bu::indexChunks then builds a local octree from each.
    distribute(m_b, fileInfos, grid, lut, progress);

    // Chunking step 3, re-chunk: the budget is soft — a single count-grid cell over the budget
    // becomes its own chunk whole, so a degenerate cell can exceed the Indexing phase's 32-bit
    // index space. Recursively split any chunk past that ceiling (packed-point read/route, no
    // decode) back down to the planning budget. Chunks merely over the budget are left alone —
    // the chunk-local build splits those in RAM. Costs nothing when nothing crosses the trigger,
    // the normal case.
    size_t split = rechunkOversized(m_b.opts.tempDir, m_b.bounds, m_b.pointSize, lut,
        cellCounts, RechunkTriggerPoints, target);
    if (split && m_b.opts.progressDebug)
        std::cerr << "[chunking] re-chunked " << split << " oversized chunk(s)\n";
}

} // namespace chunker
} // namespace untwine
