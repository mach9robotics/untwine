/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Unit tests for epf::Grid -- the octree-depth + point-to-voxel math.      *
 *   calcLevel() is the heart of the prep step (how deep the tree starts);    *
 *   key() is the heart of the EPF binning step (which voxel a point lands    *
 *   in) and is reused by reprocess at finer levels.                         *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <pdal/util/Bounds.hpp>

#include "pointio/Grid.hpp"
#include "pointio/EpfTypes.hpp"

using namespace untwine;
using namespace untwine::epf;

namespace
{
// A unit-ish cube so cell sizes are easy to reason about.
pdal::BOX3D cube(double side = 16.0)
{
    return pdal::BOX3D(0, 0, 0, side, side, side);
}
}

// ---- prep step: calcLevel() -------------------------------------------------

// A cloud at or below the per-node cap needs no subdivision: level 0.
TEST(Grid, calc_level_zero_when_under_cap)
{
    Grid g(cube(), MaxPointsPerNode, /*level*/ -1, /*cubic*/ false);
    // Auto-leveling clamps the *grid* to at least level 1 (resetLevel), but the
    // computed octree depth before clamping is 0 for a cloud at the cap.
    EXPECT_EQ(g.calcLevel(MaxPointsPerNode, false), 0);
    EXPECT_EQ(g.calcLevel(MaxPointsPerNode - 1, false), 0);
}

// Non-cubic: each level divides the count by 8 (one octree subdivision).
TEST(Grid, calc_level_noncubic_divides_by_eight)
{
    Grid g(cube(), 1, -1, false);
    // Just over the cap -> one /8 subdivision.
    EXPECT_EQ(g.calcLevel(MaxPointsPerNode + 1, false), 1);
    // count/8 still over the cap -> a second subdivision (100000*64 -> /8 ->
    // 800000 -> /8 -> 100000, stop).
    EXPECT_EQ(g.calcLevel(MaxPointsPerNode * 64, false), 2);
    // Monotonic non-decreasing in point count.
    EXPECT_GE(g.calcLevel(MaxPointsPerNode * 64, false),
              g.calcLevel(MaxPointsPerNode * 8, false));
}

// Cubic over a cube behaves identically to non-cubic (all three sides divide).
TEST(Grid, calc_level_cubic_cube_equals_noncubic)
{
    Grid g(cube(), 1, -1, true);
    uint64_t n = MaxPointsPerNode * 100;
    EXPECT_EQ(g.calcLevel(n, true), g.calcLevel(n, false));
}

// Cubic over a flat (non-cube) extent only halves along the long axis per
// level, so it needs *more* levels than the /8 non-cubic path.
TEST(Grid, calc_level_cubic_flat_extent)
{
    // Long in x, thin in y/z: only x >= side each level -> count halves, not /8.
    pdal::BOX3D flat(0, 0, 0, 1000, 0.001, 0.001);
    Grid g(flat, 1, -1, true);
    EXPECT_GT(g.calcLevel(MaxPointsPerNode * 8, true),
              g.calcLevel(MaxPointsPerNode * 8, false));
}

// ---- EPF step: resetLevel() + key() ----------------------------------------

// The grid is floored at level 1 even if asked for 0 (sampling needs >=1).
TEST(Grid, reset_level_minimum_one)
{
    Grid g(cube(), 1, 0, false);
    EXPECT_EQ(g.maxLevel(), 1);
}

// A point at the min corner lands in voxel (0,0,0); a point just inside the max
// corner lands in the last cell (gridSize-1); the level tags the grid depth.
TEST(Grid, key_corners)
{
    Grid g(cube(16.0), 1, /*level*/ 2, false);  // gridSize = 4, cell = 4.0
    EXPECT_EQ(g.maxLevel(), 2);

    EXPECT_EQ(g.key(0, 0, 0), VoxelKey(0, 0, 0, 2));
    EXPECT_EQ(g.key(15.999, 15.999, 15.999), VoxelKey(3, 3, 3, 2));
    EXPECT_EQ(g.key(5.0, 9.0, 13.0), VoxelKey(1, 2, 3, 2));
}

// Out-of-bounds coordinates clamp into the valid [0, gridSize-1] range rather
// than producing negative or oversized indices.
TEST(Grid, key_clamps_out_of_bounds)
{
    Grid g(cube(16.0), 1, 2, false);  // gridSize = 4
    EXPECT_EQ(g.key(-100, -100, -100), VoxelKey(0, 0, 0, 2));
    EXPECT_EQ(g.key(1000, 1000, 1000), VoxelKey(3, 3, 3, 2));
}

// Every key an in-bounds point can produce stays within the grid.
TEST(Grid, key_always_in_range)
{
    const int level = 4;
    const int gridSize = 1 << level;
    Grid g(cube(16.0), 1, level, false);
    for (double x = 0; x <= 16.0; x += 0.5)
        for (double y = 0; y <= 16.0; y += 4.0)
        {
            VoxelKey k = g.key(x, y, 8.0);
            EXPECT_GE(k.x(), 0);
            EXPECT_LT(k.x(), gridSize);
            EXPECT_GE(k.y(), 0);
            EXPECT_LT(k.y(), gridSize);
            EXPECT_EQ(k.level(), level);
        }
}
