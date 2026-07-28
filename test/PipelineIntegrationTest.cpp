/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Integration tests for the *seams* between pipeline steps. These exercise *
 *   the real algorithmic units (Grid, VoxelKey, Stats) wired together the    *
 *   way the orchestration code wires them, without spinning up the threaded  *
 *   FileProcessor / Reprocessor / PyramidManager machinery:                  *
 *                                                                           *
 *     prep -> EPF        : a grid sized from prep maps every in-bounds point  *
 *                          into a valid voxel; binning loses no points.       *
 *     EPF  -> reprocess  : an over-cap voxel re-splits into deeper children   *
 *                          that conserve every point.                         *
 *     reprocess -> BU    : deeper child keys roll up to the exact parent key  *
 *                          (octree consistency the BU walk depends on).       *
 *     BU   -> output      : per-node Stats merged up the tree equal the       *
 *                          single-pass whole-cloud Stats in the header.       *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <map>
#include <vector>

#include <pdal/util/Bounds.hpp>

#include "pointio/Grid.hpp"
#include "pointio/EpfTypes.hpp"
#include "untwine/VoxelKey.hpp"
#include "bu/Stats.hpp"

using namespace untwine;
using namespace untwine::epf;

namespace
{
struct Pt { double x, y, z; };

// A deterministic spread of points strictly inside [0,side)^3.
std::vector<Pt> syntheticCloud(double side, int per_axis)
{
    std::vector<Pt> pts;
    const double step = side / per_axis;
    for (int i = 0; i < per_axis; ++i)
        for (int j = 0; j < per_axis; ++j)
            for (int k = 0; k < per_axis; ++k)
                pts.push_back({ (i + 0.5) * step, (j + 0.5) * step, (k + 0.5) * step });
    return pts;
}

// Ascend `levels` parents from a key.
VoxelKey ancestor(VoxelKey k, int levels)
{
    for (int i = 0; i < levels; ++i)
        k = k.parent();
    return k;
}
}

// ---- prep -> EPF ------------------------------------------------------------
// The depth prep picks (calcLevel) and the grid EPF bins with (key) compose:
// every point inside the prep bounds bins into an in-range voxel at the grid's
// level, and binning is conservative (no point dropped, none invented).
TEST(Integration, prep_to_epf_binning_conserves_points)
{
    const double side = 64.0;
    pdal::BOX3D bounds(0, 0, 0, side, side, side);
    std::vector<Pt> cloud = syntheticCloud(side, 16);  // 4096 points

    // prep: choose the octree level from the point count.
    Grid grid(bounds, cloud.size(), /*level*/ -1, /*cubic*/ false);
    const int gridSize = 1 << grid.maxLevel();

    // EPF: bin every point.
    std::map<VoxelKey, size_t> bins;
    for (const Pt& p : cloud)
    {
        VoxelKey k = grid.key(p.x, p.y, p.z);
        EXPECT_EQ(k.level(), grid.maxLevel());
        EXPECT_GE(k.x(), 0); EXPECT_LT(k.x(), gridSize);
        EXPECT_GE(k.y(), 0); EXPECT_LT(k.y(), gridSize);
        EXPECT_GE(k.z(), 0); EXPECT_LT(k.z(), gridSize);
        bins[k]++;
    }

    size_t total = 0;
    for (const auto& kv : bins)
        total += kv.second;
    EXPECT_EQ(total, cloud.size());  // conservation across the prep->EPF seam
}

// ---- EPF -> reprocess -------------------------------------------------------
// When a coarse voxel holds more than MaxPointsPerNode points, reprocess
// re-keys those same points on a finer grid (parent level + extra levels). The
// children must (a) all descend from the original voxel and (b) conserve every
// point -- the same guarantee Reprocessor::run() relies on.
TEST(Integration, epf_to_reprocess_resplit_conserves_and_descends)
{
    const double side = 8.0;
    pdal::BOX3D bounds(0, 0, 0, side, side, side);

    // Coarse EPF grid at level 1: cell size 4, so points in [0,4)^3 all land in
    // the single voxel (0,0,0,1). Cluster the cloud there to force a re-split.
    Grid coarse(bounds, 1, /*level*/ 1, false);
    std::vector<Pt> cluster;
    for (int i = 0; i < 32; ++i)
        for (int j = 0; j < 32; ++j)
            cluster.push_back({ i * (3.9 / 32), j * (3.9 / 32), 1.0 });

    VoxelKey parent = coarse.key(cluster.front().x, cluster.front().y, cluster.front().z);
    ASSERT_EQ(parent, VoxelKey(0, 0, 0, 1));
    for (const Pt& p : cluster)
        ASSERT_EQ(coarse.key(p.x, p.y, p.z), parent);  // all in one coarse voxel

    // reprocess: split with `extraLevels` extra subdivisions (as Reprocessor
    // does via resetLevel(parentLevel + levels)).
    const int extraLevels = 2;
    Grid fine = coarse;
    fine.resetLevel(parent.level() + extraLevels);

    std::map<VoxelKey, size_t> children;
    for (const Pt& p : cluster)
    {
        VoxelKey c = fine.key(p.x, p.y, p.z);
        EXPECT_EQ(c.level(), parent.level() + extraLevels);
        EXPECT_EQ(ancestor(c, extraLevels), parent);   // descends from the parent voxel
        children[c]++;
    }

    size_t total = 0;
    for (const auto& kv : children)
        total += kv.second;
    EXPECT_EQ(total, cluster.size());                   // conservation across the re-split
    EXPECT_GT(children.size(), 1u);                     // actually subdivided
}

// ---- reprocess -> BU --------------------------------------------------------
// BU builds the tree leaf-to-root by repeatedly taking parent(). For that to
// land points correctly, a key produced at a deep level must, after ascending
// the level delta, equal the key the same point gets at the shallow level.
// This is what links reprocess's deeper keys back to their EPF ancestors.
TEST(Integration, reprocess_to_bu_octree_consistency)
{
    const double side = 100.0;
    pdal::BOX3D bounds(0, 0, 0, side, side, side);
    std::vector<Pt> cloud = syntheticCloud(side, 12);

    const int shallow = 3;
    const int deep = 6;
    Grid shallowGrid(bounds, 1, shallow, false);
    Grid deepGrid(bounds, 1, deep, false);

    for (const Pt& p : cloud)
    {
        VoxelKey shallowKey = shallowGrid.key(p.x, p.y, p.z);
        VoxelKey deepKey = deepGrid.key(p.x, p.y, p.z);
        // Ascend the deep key back up to the shallow level == EPF/BU rollup.
        EXPECT_EQ(ancestor(deepKey, deep - shallow), shallowKey);
        // And all the way to the root.
        EXPECT_EQ(ancestor(deepKey, deep), VoxelKey(0, 0, 0, 0));
    }
}

// ---- BU -> output -----------------------------------------------------------
// The BU step accumulates Stats per node and merges them up the tree; the root
// merge is what gets written to the COPC header. Merging per-node partial stats
// must equal a single pass over the whole cloud.
TEST(Integration, bu_to_output_stats_rollup)
{
    // Partition values into "nodes" exactly as the tree would.
    std::vector<std::vector<double>> nodes {
        { 1, 2, 3, 4 },
        { 10, 12, 14 },
        { -5, -5, 0, 5, 5 },
        { 100 },
    };

    Stats whole("Z", Stats::NoEnum);
    for (const auto& node : nodes)
        for (double v : node)
            whole.insert(v);

    // Per-node stats, then merge up (order-independent fold).
    Stats rolled("Z", Stats::NoEnum);
    for (const auto& node : nodes)
    {
        Stats nodeStats("Z", Stats::NoEnum);
        for (double v : node)
            nodeStats.insert(v);
        ASSERT_TRUE(rolled.merge(nodeStats));
    }

    EXPECT_EQ(rolled.count(), whole.count());
    EXPECT_DOUBLE_EQ(rolled.minimum(), whole.minimum());
    EXPECT_DOUBLE_EQ(rolled.maximum(), whole.maximum());
    EXPECT_NEAR(rolled.average(), whole.average(), 1e-9);
    EXPECT_NEAR(rolled.populationStddev(), whole.populationStddev(), 1e-9);
}
