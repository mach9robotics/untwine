/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Unit tests for the chunk-local BU build (bu/ChunkPlan): the adaptive     *
 *   chunk-root planner and the top-down octant partition that buildHierarchy *
 *   uses. These are the pure-logic seams of the chunked build; the threaded  *
 *   ChunkBuilder/PyramidManager orchestration is exercised end-to-end by the *
 *   binary test.                                                            *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <unordered_map>
#include <vector>

#include <pdal/util/Bounds.hpp>

#include "bu/ChunkPlan.hpp"
#include "untwine/VoxelKey.hpp"

using namespace untwine;
using namespace untwine::bu;

namespace
{

// Is `node` a descendant-or-self of `root`?
bool descendsFrom(VoxelKey node, VoxelKey root)
{
    if (node.level() < root.level())
        return false;
    while (node.level() > root.level())
        node = node.parent();
    return node == root;
}

// Total leaf points in `root`'s subtree.
uint64_t subtreeCount(const std::unordered_map<VoxelKey, uint64_t>& leaves, VoxelKey root)
{
    uint64_t t = 0;
    for (const auto& p : leaves)
        if (descendsFrom(p.first, root))
            t += p.second;
    return t;
}

// nPerAxis^3 leaves at `level`, each with `count` points.
std::unordered_map<VoxelKey, uint64_t> makeLeaves(int level, int nPerAxis, uint64_t count)
{
    std::unordered_map<VoxelKey, uint64_t> m;
    for (int x = 0; x < nPerAxis; ++x)
        for (int y = 0; y < nPerAxis; ++y)
            for (int z = 0; z < nPerAxis; ++z)
                m[VoxelKey(x, y, z, level)] = count;
    return m;
}

// Add an n^3 block of leaves at `level` with its corner at (x0,y0,z0), each with `count` points.
// Used to place dense leaf clusters at arbitrary offsets (far apart) in a sparse bbox.
void addBlock(std::unordered_map<VoxelKey, uint64_t>& m, int level,
    int x0, int y0, int z0, int n, uint64_t count)
{
    for (int x = 0; x < n; ++x)
        for (int y = 0; y < n; ++y)
            for (int z = 0; z < n; ++z)
                m[VoxelKey(x0 + x, y0 + y, z0 + z, level)] = count;
}

// How many of `roots` is `leaf` a descendant of?
int rootAncestors(const std::vector<VoxelKey>& roots, VoxelKey leaf)
{
    int n = 0;
    for (const VoxelKey& r : roots)
        if (descendsFrom(leaf, r))
            ++n;
    return n;
}

} // namespace

// The planner must carve the leaves into a valid partition: chunk roots that are disjoint,
// each within the point budget, below the global root, and never double-covering a leaf.
TEST(ChunkPlan, roots_are_a_valid_partition)
{
    auto leaves = makeLeaves(/*level=*/4, /*nPerAxis=*/4, /*count=*/1000);  // 64 leaves, 64k pts
    const uint64_t target = 8000;
    std::vector<VoxelKey> roots = planChunkRoots(leaves, target);
    ASSERT_FALSE(roots.empty());

    // (1) Chunk roots are never the global root (always level >= 1), so the cap has work.
    for (const VoxelKey& r : roots)
        EXPECT_GE(r.level(), 1);

    // (2) Disjoint: no root is an ancestor/descendant of another.
    for (size_t i = 0; i < roots.size(); ++i)
        for (size_t j = 0; j < roots.size(); ++j)
            if (i != j)
                EXPECT_FALSE(descendsFrom(roots[i], roots[j]))
                    << roots[i] << " nested under " << roots[j];

    // (3) Each chunk root's subtree fits the budget.
    for (const VoxelKey& r : roots)
        EXPECT_LE(subtreeCount(leaves, r), target);

    // (4) No leaf is claimed by more than one internal root.
    for (const auto& p : leaves)
        EXPECT_LE(rootAncestors(roots, p.first), 1) << "leaf " << p.first;
}

// The budget controls how deep the chunk-root cut lands. (Root count is NOT monotonic in the
// budget: a too-small budget descends past internal nodes onto the leaves, which then become
// standalone single-leaf roots rather than internal ones.) Here 64 leaves all sit under one
// level-1 node, so a per-2x2x2-block budget cuts at level 3 (8 roots) while a whole-cloud
// budget cuts at level 1 (1 root) -- a larger budget => shallower, fewer roots.
TEST(ChunkPlan, budget_controls_cut_level)
{
    auto leaves = makeLeaves(/*level=*/4, /*nPerAxis=*/4, /*count=*/1000);  // 64 leaves, 64k pts

    std::vector<VoxelKey> deep = planChunkRoots(leaves, 8000);   // one 2x2x2 block per chunk
    EXPECT_EQ(deep.size(), 8u);
    for (const VoxelKey& r : deep)
        EXPECT_EQ(r.level(), 3);

    std::vector<VoxelKey> shallow = planChunkRoots(leaves, 64000);  // whole cloud in one chunk
    EXPECT_EQ(shallow.size(), 1u);
    EXPECT_EQ(shallow.front().level(), 1);
}

// With the on-disk reprocessor retired, an oversized leaf can't be pre-split, so when no
// ancestor fits the budget the leaf itself becomes a buildable chunk root (buildHierarchy
// splits it in RAM). A leaf that DOES fit its ancestors is rolled up into a shallow internal
// root instead -- a lone leaf is never left standalone here, only leaves in an over-budget
// neighborhood are (covered by roots_are_a_valid_partition).
TEST(ChunkPlan, oversized_leaf_is_buildable)
{
    VoxelKey leaf(1, 2, 3, 5);

    std::unordered_map<VoxelKey, uint64_t> big{ { leaf, 1000000 } };
    std::vector<VoxelKey> roots = planChunkRoots(big, 100);
    ASSERT_EQ(roots.size(), 1u);
    EXPECT_EQ(roots.front(), leaf);

    std::unordered_map<VoxelKey, uint64_t> small{ { leaf, 50 } };
    std::vector<VoxelKey> sroots = planChunkRoots(small, 100);
    ASSERT_EQ(sroots.size(), 1u);
    EXPECT_EQ(sroots.front().level(), 1);  // rolled up to a shallow internal root, not the leaf
}

// When the whole cloud fits one budget, the root itself would qualify but is excluded, so the
// cut still lands at level >= 1 and every leaf ends up under exactly one root.
TEST(ChunkPlan, whole_cloud_splits_below_root)
{
    auto leaves = makeLeaves(/*level=*/3, /*nPerAxis=*/2, /*count=*/10);  // 8 leaves under (0,0,0,1)
    std::vector<VoxelKey> roots = planChunkRoots(leaves, 1000000);
    ASSERT_FALSE(roots.empty());
    for (const VoxelKey& r : roots)
        EXPECT_GE(r.level(), 1);
    for (const auto& p : leaves)
        EXPECT_EQ(rootAncestors(roots, p.first), 1) << "leaf " << p.first;
}

// childOctant indexes the 8 children with the same bit layout buildHierarchy relies on:
// bit0->x, bit1->y, bit2->z, split at the node midpoint, midpoint going to the upper side.
TEST(ChunkBuild, child_octant_bit_layout)
{
    pdal::BOX3D b(0, 0, 0, 8, 8, 8);  // midpoint (4,4,4)
    EXPECT_EQ(childOctant(1, 1, 1, b), 0);  // all low
    EXPECT_EQ(childOctant(6, 1, 1, b), 1);  // x high -> bit0
    EXPECT_EQ(childOctant(1, 6, 1, b), 2);  // y high -> bit1
    EXPECT_EQ(childOctant(1, 1, 6, b), 4);  // z high -> bit2
    EXPECT_EQ(childOctant(6, 6, 6, b), 7);  // all high
    EXPECT_EQ(childOctant(4, 4, 4, b), 7);  // exactly on midpoint -> upper side
}

// The octant index must agree with VoxelKey::child(): a point in octant `dir` belongs to
// child(dir), whose low coordinate bits encode the same x/y/z halves. A mismatch here would
// route points to the wrong child during the in-chunk build.
TEST(ChunkBuild, child_octant_matches_voxelkey_child)
{
    pdal::BOX3D b(0, 0, 0, 8, 8, 8);
    VoxelKey parent(0, 0, 0, 0);
    for (int dir = 0; dir < 8; ++dir)
    {
        double x = (dir & 1) ? 6.0 : 2.0;
        double y = (dir & 2) ? 6.0 : 2.0;
        double z = (dir & 4) ? 6.0 : 2.0;
        EXPECT_EQ(childOctant(x, y, z, b), dir);

        VoxelKey c = parent.child(dir);
        EXPECT_EQ(c.x() & 1, (dir >> 0) & 1);
        EXPECT_EQ(c.y() & 1, (dir >> 1) & 1);
        EXPECT_EQ(c.z() & 1, (dir >> 2) & 1);
    }
}

// Models two_traj AFTER the on-disk reprocessor refines its coarse voxels: two dense clusters of
// fine leaves sit far apart with a large empty gap between them (a huge, mostly-empty bbox). This
// is the case the rebinning fix restores. The clusters diverge immediately below the root -- one in
// octant 0, the other in octant 7 of level 1 -- so no shallow node spans both, and the empty gap
// can't collapse them into one chunk. The planner must carve many balanced chunks per cluster,
// which is what gives BU its parallelism.
TEST(ChunkPlan, sparse_far_apart_clusters_split_into_many_balanced_chunks)
{
    std::unordered_map<VoxelKey, uint64_t> leaves;
    // Two 4x4x4 leaf blocks at level 6 (coords 0..63), at opposite corners of the bbox.
    addBlock(leaves, /*level=*/6, /*x0,y0,z0=*/0, 0, 0,    /*n=*/4, /*count=*/1000);  // 64k pts
    addBlock(leaves, /*level=*/6, /*x0,y0,z0=*/60, 60, 60, /*n=*/4, /*count=*/1000);  // 64k pts

    const uint64_t target = 8000;
    std::vector<VoxelKey> roots = planChunkRoots(leaves, target);

    // Many chunks, not ~3: each 2x2x2 sub-block (8 leaves = 8000 pts) becomes one chunk root, so
    // each cluster yields 8 -> 16 total. The point vs the regression is "far more than the handful
    // of coarse occupied cells", not the exact 16.
    EXPECT_GT(roots.size(), 2u);
    EXPECT_EQ(roots.size(), 16u);

    // Balanced: every chunk fits the budget (none is an oversized coarse voxel).
    for (const VoxelKey& r : roots)
        EXPECT_LE(subtreeCount(leaves, r), target) << "oversized chunk " << r;

    // Valid disjoint partition: every leaf is under exactly one chunk root.
    for (const auto& p : leaves)
        EXPECT_EQ(rootAncestors(roots, p.first), 1) << "leaf " << p.first;

    // Both far-apart clusters are independently subdivided, not collapsed by the empty gap.
    int clusterA = 0, clusterB = 0;
    for (const VoxelKey& r : roots)
    {
        if (descendsFrom(r, VoxelKey(0, 0, 0, 1)))
            ++clusterA;
        if (descendsFrom(r, VoxelKey(1, 1, 1, 1)))
            ++clusterB;
    }
    EXPECT_EQ(clusterA, 8);
    EXPECT_EQ(clusterB, 8);
}

// The pre-reprocess state that produced the 244 s regression: the same two far-apart regions, but
// each is a single coarse over-budget leaf (what EPF's coarse grid yields before reprocess refines
// it). planChunkRoots can only return one buildable root per occupied voxel, so BU runs with as
// many independent chunks as there are coarse cells -- here just 2 -> ~2-way parallel. This is the
// starvation the rebinning fix removes; contrast the many balanced chunks above.
TEST(ChunkPlan, coarse_far_apart_leaves_starve_parallelism)
{
    std::unordered_map<VoxelKey, uint64_t> leaves{
        { VoxelKey(0, 0, 0, 2), 64000 },
        { VoxelKey(3, 3, 3, 2), 64000 },
    };
    const uint64_t target = 8000;
    std::vector<VoxelKey> roots = planChunkRoots(leaves, target);

    // One oversized chunk per occupied coarse voxel -> minimal parallelism, unbalanced chunks.
    EXPECT_EQ(roots.size(), 2u);
    for (const VoxelKey& r : roots)
        EXPECT_GT(subtreeCount(leaves, r), target) << "chunk " << r << " should be over budget";
}

// localGridDepth: per-chunk local count-grid depth. 8^d should cover ceil(N/budget) with one
// level of margin, capped at maxDepthBelowRoot.
TEST(LocalGridDepth, tiers_and_cap)
{
    const uint64_t budget = 100000;
    // N <= budget -> ratio 1 -> base depth 1 + 1 margin = 2.
    EXPECT_EQ(localGridDepth(budget, budget, 10), 2);
    EXPECT_EQ(localGridDepth(1, budget, 10), 2);
    // ratio 50 -> 8^2=64 covers -> depth 2 + margin = 3.
    EXPECT_EQ(localGridDepth(50 * budget, budget, 10), 3);
    // ratio 1000 -> 8^4=4096 covers (8^3=512 < 1000) -> depth 4 + margin = 5.
    EXPECT_EQ(localGridDepth(1000 * budget, budget, 10), 5);
    // Cap at maxDepthBelowRoot.
    EXPECT_EQ(localGridDepth(50 * budget, budget, 2), 2);
    EXPECT_GE(localGridDepth(50 * budget, budget, 1), 1);
    EXPECT_LE(localGridDepth(1000 * budget, budget, 3), 3);
}
