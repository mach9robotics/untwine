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

// Trigger for the re-chunk pass: a chunk over this many points is split after distribute. Set to
// the hard ceiling of the chunk-local build's 32-bit index space, so re-chunking is purely a
// correctness repair for cells the planner couldn't resolve — chunks between the planning budget
// and this trigger are left alone (the chunk-local build handles them fine given RAM; splitting
// them costs a full rewrite of their temp bytes). Once triggered, a chunk is split down to the
// planning budget, restoring normal per-chunk RAM and parallelism.
constexpr uint64_t RechunkTriggerPoints = 2'147'483'647;

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

// Fallback worker-thread count when std::thread::hardware_concurrency() reports 0.
constexpr unsigned DefaultWorkerThreads = 8;

// Points per PDAL streaming batch (FixedPointTable capacity) during the count and distribute
// decodes. Matches the value EPF's FileProcessor uses.
constexpr std::size_t StreamPointBatch = 1000;

} // namespace chunker
} // namespace untwine
