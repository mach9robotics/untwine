/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Pure-logic unit tests for the chunker (counting-sort) front-end:         *
 *   chunker::countGridLevel (count-grid tier selection) and                  *
 *   chunker::buildChunkLut, which partitions the count-grid histogram into     *
 *   chunks: the cell->chunk map the distribute pass routes points through.      *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <unordered_map>
#include <unordered_set>

#include "chunker/ChunkLut.hpp"
#include "chunker/CountPass.hpp"
#include "untwine/VoxelKey.hpp"

using namespace untwine;
using namespace untwine::chunker;

namespace
{

bool descendsFrom(VoxelKey node, VoxelKey root)
{
    if (node.level() < root.level())
        return false;
    while (node.level() > root.level())
        node = node.parent();
    return node == root;
}

std::unordered_map<VoxelKey, uint64_t> makeCells(int level, int nPerAxis, uint64_t count)
{
    std::unordered_map<VoxelKey, uint64_t> m;
    for (int x = 0; x < nPerAxis; ++x)
        for (int y = 0; y < nPerAxis; ++y)
            for (int z = 0; z < nPerAxis; ++z)
                m[VoxelKey(x, y, z, level)] = count;
    return m;
}

} // namespace

// Every occupied count-cell maps to exactly one chunk; that chunk is an ancestor-or-self of the
// cell and is in the root set; and the roots are disjoint.
TEST(ChunkLut, every_cell_maps_to_one_covering_chunk)
{
    auto cells = makeCells(/*level=*/4, /*nPerAxis=*/4, /*count=*/1000);  // 64 cells, 64k pts
    ChunkLut lut = buildChunkLut(cells, 8000);
    std::unordered_set<VoxelKey> rootSet(lut.roots.begin(), lut.roots.end());

    EXPECT_EQ(lut.cellToRoot.size(), cells.size());
    for (const auto& c : cells)
    {
        auto it = lut.cellToRoot.find(c.first);
        ASSERT_NE(it, lut.cellToRoot.end()) << "unmapped cell " << c.first;
        EXPECT_TRUE(descendsFrom(c.first, it->second))
            << "cell " << c.first << " not under chunk " << it->second;
        EXPECT_TRUE(rootSet.count(it->second)) << "chunk " << it->second << " missing from roots";
    }

    for (size_t i = 0; i < lut.roots.size(); ++i)
        for (size_t j = 0; j < lut.roots.size(); ++j)
            if (i != j)
                EXPECT_FALSE(descendsFrom(lut.roots[i], lut.roots[j]))
                    << lut.roots[i] << " nested under " << lut.roots[j];
}

// Two sibling cells whose combined count exceeds the budget can't merge into their parent, so
// each becomes its own chunk at the count-grid level.
TEST(ChunkLut, unmergeable_cells_become_their_own_chunks)
{
    std::unordered_map<VoxelKey, uint64_t> cells{
        { VoxelKey(0, 0, 0, 5), 6000 },
        { VoxelKey(1, 0, 0, 5), 6000 },
    };
    ChunkLut lut = buildChunkLut(cells, 8000);  // 6000 + 6000 = 12000 > 8000
    EXPECT_EQ(lut.cellToRoot.at(VoxelKey(0, 0, 0, 5)), VoxelKey(0, 0, 0, 5));
    EXPECT_EQ(lut.cellToRoot.at(VoxelKey(1, 0, 0, 5)), VoxelKey(1, 0, 0, 5));
}

// A target large enough to hold the whole occupied region merges every cell up into a single
// chunk root (the shallowest non-root ancestor whose subtree fits). The 4x4x4 block at level 4
// all descends from the level-1 node (0,0,0,1).
TEST(ChunkLut, large_target_merges_cells_to_one_root)
{
    auto cells = makeCells(/*level=*/4, /*nPerAxis=*/4, /*count=*/1000);  // 64 cells, 64k pts
    ChunkLut lut = buildChunkLut(cells, 1000000);                         // target >> total
    ASSERT_EQ(lut.roots.size(), 1u);
    EXPECT_EQ(lut.roots[0], VoxelKey(0, 0, 0, 1));
    for (const auto& c : cells)
        EXPECT_EQ(lut.cellToRoot.at(c.first), VoxelKey(0, 0, 0, 1));
}

// A single cell whose own count already exceeds the target can't merge anywhere; it is a buildable
// chunk root in its own right (ChunkBuilder later splits it in RAM).
TEST(ChunkLut, oversized_single_cell_is_its_own_root)
{
    std::unordered_map<VoxelKey, uint64_t> cells{ { VoxelKey(5, 5, 5, 4), 200000 } };
    ChunkLut lut = buildChunkLut(cells, 8000);  // 200000 > 8000
    ASSERT_EQ(lut.roots.size(), 1u);
    EXPECT_EQ(lut.roots[0], VoxelKey(5, 5, 5, 4));
    EXPECT_EQ(lut.cellToRoot.at(VoxelKey(5, 5, 5, 4)), VoxelKey(5, 5, 5, 4));
}

// Count-grid tier boundaries (128^3 / 256^3 / 512^3 = levels 7 / 8 / 9).
TEST(CountGrid, level_tier_boundaries)
{
    EXPECT_EQ(countGridLevel(0), 7);
    EXPECT_EQ(countGridLevel(99999999ULL), 7);
    EXPECT_EQ(countGridLevel(100000000ULL), 8);   // boundary is exclusive on the low side
    EXPECT_EQ(countGridLevel(499999999ULL), 8);
    EXPECT_EQ(countGridLevel(500000000ULL), 9);
    EXPECT_EQ(countGridLevel(5000000000ULL), 9);
}
