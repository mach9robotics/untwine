/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                      *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 ****************************************************************************/

#pragma once

#include <cstddef>
#include <cstdint>

// Tunable constants for the chunker (counting-sort) pipeline, collected so the policy values are
// named and live in one place rather than as literals at their use sites.

namespace untwine
{
namespace chunker
{

// Default per-chunk point budget used when --max_chunk_points is 0/unset. This is the PLANNING
// threshold: chunk boundaries merge up the octree until one level higher would exceed it. Chosen
// from profiling: large enough to keep the merge phase small, small enough for good per-chunk
// parallelism (both extremes lose — tiny budgets bloat the merge, huge ones starve distribute and
// the indexing pool; 2^31-1 as budget cost +30% wall on two_traj).
constexpr uint64_t DefaultMaxChunkPoints = 5'000'000;

// Trigger for the re-chunk pass: a chunk over this many points is split after distribute. A pure
// performance bound, not a correctness one (the chunk-local build picks a 64-bit index width for
// over-2^31 chunks and handles any size, bounded by RAM): past the trigger, restoring per-chunk
// RAM and parallelism outweighs rewriting the chunk's temp bytes. Chunks between the planning
// budget and this trigger are left to the in-RAM build. Once triggered, a chunk is split down to
// the planning budget. Chosen from a 10/20/50/100M sweep on medium_mx9 (one 382M-point count
// cell, 16 GB box): 10-50M tie at the optimum (-46% wall vs no split), 100M forfeits part of the
// win through the merge. Overridable via UNTWINE_RECHUNK_TRIGGER.
constexpr uint64_t RechunkTriggerPoints = 20'000'000;

// Count-grid resolution tiers (see countGridLevel). The count pass histograms points into an
// octree grid this many levels deep, finer for larger clouds so cells stay small relative to the
// chunk budget. A cloud with fewer than ...SmallMaxPoints uses the small grid, fewer than
// ...MediumMaxPoints the medium grid, otherwise the large grid.
constexpr uint64_t CountGridSmallMaxPoints = 100'000'000;
constexpr uint64_t CountGridMediumMaxPoints = 500'000'000;
constexpr int CountGridLevelSmall = 7;   // 128^3
constexpr int CountGridLevelMedium = 8;  // 256^3
constexpr int CountGridLevelLarge = 9;   // 512^3

// Depth cap for the re-chunk pass's recursive split of oversized chunks. A coincident point
// cluster occupies one octant at every level, so without a floor the split would recurse forever;
// matches the MaxBuildLevel guard the Indexing phase uses for the same pathology.
constexpr int MaxRechunkLevel = 30;

// Bounds for one re-chunk split step. MaxSplitStepDepth caps how many levels a single
// read+write pass descends (8^d output files and staging buffers open at once: depth 3 = up to
// 512). SplitStageBytes is the per-output staging buffer, batching records into block-sized
// writes; worst case memory per split task = 8^MaxSplitStepDepth * SplitStageBytes = 512 MB.
constexpr int MaxSplitStepDepth = 3;
constexpr std::size_t SplitStageBytes = 1024 * 1024;

// Fallback worker-thread count when std::thread::hardware_concurrency() reports 0.
constexpr unsigned DefaultWorkerThreads = 8;

// Points per PDAL streaming batch (FixedPointTable capacity) during the count and distribute
// decodes. Matches the value EPF's FileProcessor uses.
constexpr std::size_t StreamPointBatch = 1000;

} // namespace chunker
} // namespace untwine
