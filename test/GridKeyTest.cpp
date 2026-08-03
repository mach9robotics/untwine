/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Unit tests for GridKey -- the (i,j,k) cell index used by the BU          *
 *   Processor's sampling grid (one accepted point per occupied cell).       *
 ****************************************************************************/

#include <gtest/gtest.h>

#include "untwine/GridKey.hpp"

using namespace untwine;

// i/j/k pack into a single int as (i<<16)|(j<<8)|k and unpack losslessly for
// values that fit in a byte each (the grid is CellCount=128 per axis).
TEST(GridKey, ijk_round_trip)
{
    GridKey k(127, 63, 1);
    EXPECT_EQ(k.i(), 127);
    EXPECT_EQ(k.j(), 63);
    EXPECT_EQ(k.k(), 1);

    GridKey zero(0, 0, 0);
    EXPECT_EQ(zero.i(), 0);
    EXPECT_EQ(zero.j(), 0);
    EXPECT_EQ(zero.k(), 0);
}

// The packed key value is the byte layout we expect.
TEST(GridKey, packed_key_layout)
{
    GridKey k(0x12, 0x34, 0x56);
    EXPECT_EQ(k.key(), (0x12 << 16) | (0x34 << 8) | 0x56);
}

TEST(GridKey, equality)
{
    EXPECT_EQ(GridKey(5, 6, 7), GridKey(5, 6, 7));
    EXPECT_FALSE(GridKey(5, 6, 7) == GridKey(5, 6, 8));
}

// The hash is just the hash of the packed key, so equal keys hash equally.
TEST(GridKey, hash_matches_key)
{
    GridKey k(10, 20, 30);
    EXPECT_EQ(std::hash<GridKey>()(k), std::hash<int>()(k.key()));
}
