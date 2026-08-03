/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Unit tests for VoxelKey -- the octree node key that is the common        *
 *   currency of every pipeline step: EPF bins points into VoxelKeys,        *
 *   reprocess splits a VoxelKey into deeper children, and BU walks the       *
 *   tree from leaves to root via parent()/child().                          *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <unordered_set>

#include "untwine/VoxelKey.hpp"

using namespace untwine;

// A default key is the root: origin at level 0.
TEST(VoxelKey, default_is_root)
{
    VoxelKey k;
    EXPECT_EQ(k.x(), 0);
    EXPECT_EQ(k.y(), 0);
    EXPECT_EQ(k.z(), 0);
    EXPECT_EQ(k.level(), 0);
}

// child(dir) encodes the three octant bits as: bit0->x, bit1->y, bit2->z, and
// always descends exactly one level. This is the bit layout reprocess relies on.
TEST(VoxelKey, child_bit_encoding)
{
    VoxelKey root;  // (0,0,0,0)

    // Each of the 8 children doubles the parent coordinate and ORs in the bit.
    for (int dir = 0; dir < 8; ++dir)
    {
        VoxelKey c = root.child(dir);
        EXPECT_EQ(c.level(), 1);
        EXPECT_EQ(c.x(), (dir >> 0) & 0x1);
        EXPECT_EQ(c.y(), (dir >> 1) & 0x1);
        EXPECT_EQ(c.z(), (dir >> 2) & 0x1);
    }

    // Non-origin parent: coordinates shift left by one before the bit is set.
    VoxelKey p(3, 5, 6, 4);
    VoxelKey c = p.child(0b101);  // x bit + z bit set, y bit clear
    EXPECT_EQ(c.x(), (3 << 1) | 1);
    EXPECT_EQ(c.y(), (5 << 1) | 0);
    EXPECT_EQ(c.z(), (6 << 1) | 1);
    EXPECT_EQ(c.level(), 5);
}

// parent() is the exact inverse of child() for every octant -- the invariant
// that makes the BU bottom-up rollup land each child on the right parent.
TEST(VoxelKey, parent_inverts_child)
{
    VoxelKey nodes[] = { VoxelKey(), VoxelKey(1, 0, 0, 1), VoxelKey(7, 3, 12, 5) };
    for (const VoxelKey& n : nodes)
        for (int dir = 0; dir < 8; ++dir)
            EXPECT_EQ(n.child(dir).parent(), n) << "dir=" << dir << " node=" << n;
}

// parent() of the root clamps the level at 0 rather than going negative.
TEST(VoxelKey, parent_of_root_clamps_level)
{
    VoxelKey root;
    VoxelKey p = root.parent();
    EXPECT_EQ(p.level(), 0);
    EXPECT_EQ(p, root);
}

TEST(VoxelKey, equality_and_inequality)
{
    EXPECT_EQ(VoxelKey(1, 2, 3, 4), VoxelKey(1, 2, 3, 4));
    EXPECT_NE(VoxelKey(1, 2, 3, 4), VoxelKey(1, 2, 3, 5));
    EXPECT_NE(VoxelKey(1, 2, 3, 4), VoxelKey(0, 2, 3, 4));
}

// operator< orders by x, then y, then z, then level -- used wherever keys go
// into ordered containers.
TEST(VoxelKey, ordering)
{
    EXPECT_LT(VoxelKey(0, 0, 0, 9), VoxelKey(1, 0, 0, 0));   // x dominates
    EXPECT_LT(VoxelKey(1, 0, 0, 0), VoxelKey(1, 1, 0, 0));   // then y
    EXPECT_LT(VoxelKey(1, 1, 0, 0), VoxelKey(1, 1, 1, 0));   // then z
    EXPECT_LT(VoxelKey(1, 1, 1, 0), VoxelKey(1, 1, 1, 1));   // then level
}

// The .bin temp-file names are derived from this string ("level-x-y-z"); a
// regression here would scatter points into mis-named files.
TEST(VoxelKey, to_string_format)
{
    EXPECT_EQ(VoxelKey(2, 3, 4, 1).toString(), "1-2-3-4");
    EXPECT_EQ(VoxelKey().toString(), "0-0-0-0");
}

// Distinct keys (within the documented 16-bit-per-axis assumption of the hash)
// must produce distinct hashes so the EPF/BU hash maps don't collide spuriously.
TEST(VoxelKey, hash_distinct_for_distinct_keys)
{
    std::hash<VoxelKey> h;
    std::unordered_set<size_t> seen;
    int collisions = 0;
    for (int x = 0; x < 8; ++x)
        for (int y = 0; y < 8; ++y)
            for (int z = 0; z < 8; ++z)
                for (int l = 0; l < 4; ++l)
                    if (!seen.insert(h(VoxelKey(x, y, z, l))).second)
                        collisions++;
    EXPECT_EQ(collisions, 0);
}
