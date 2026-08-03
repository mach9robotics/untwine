/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Unit tests for the dimension-classification helpers used by the prep     *
 *   step to decide which input dimensions are standard LAS fields, which     *
 *   are "extra bytes", and which collapse into the packed UntwineBits byte,  *
 *   plus the FileDimInfo serialization used to pass dim layout downstream.   *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <sstream>

#include "untwine/Common.hpp"
#include "untwine/FileDimInfo.hpp"

using namespace untwine;

// Standard LAS dimensions are not extra bytes; an unknown named dim is.
TEST(Common, is_extra_dim)
{
    EXPECT_FALSE(isExtraDim("X"));
    EXPECT_FALSE(isExtraDim("Intensity"));
    EXPECT_FALSE(isExtraDim("Red"));
    EXPECT_FALSE(isExtraDim("GpsTime"));
    EXPECT_FALSE(isExtraDim(UntwineBitsDimName));   // our synthetic packed byte

    EXPECT_TRUE(isExtraDim("Amplitude"));           // vendor extra byte
    EXPECT_TRUE(isExtraDim("MyCustomField"));
}

// The exploded class-flag / scan-bit dimensions map into the packed byte.
TEST(Common, is_untwine_bits_dim)
{
    EXPECT_TRUE(isUntwineBitsDim("Synthetic"));
    EXPECT_TRUE(isUntwineBitsDim("EdgeOfFlightLine"));
    EXPECT_TRUE(isUntwineBitsDim("ScanChannel"));
    EXPECT_FALSE(isUntwineBitsDim("X"));
    EXPECT_FALSE(isUntwineBitsDim("Intensity"));
}

// Bit positions used when packing the UntwineBits byte; unknown -> -1.
TEST(Common, untwine_bit_positions)
{
    EXPECT_EQ(getUntwineBitPos("Synthetic"), 0);
    EXPECT_EQ(getUntwineBitPos("KeyPoint"), 1);
    EXPECT_EQ(getUntwineBitPos("Withheld"), 2);
    EXPECT_EQ(getUntwineBitPos("Overlap"), 3);
    EXPECT_EQ(getUntwineBitPos("ScanChannel"), 5);
    EXPECT_EQ(getUntwineBitPos("ScanDirectionFlag"), 6);
    EXPECT_EQ(getUntwineBitPos("EdgeOfFlightLine"), 7);
    EXPECT_EQ(getUntwineBitPos("NotABitDim"), -1);
}

// FileDimInfo round-trips through its stream operators (name, type, offset),
// which is how the dim layout is serialized for the downstream phases.
TEST(Common, file_dim_info_round_trip)
{
    FileDimInfo in("Intensity", pdal::Dimension::Type::Unsigned16);
    in.offset = 12;

    std::stringstream ss;
    ss << in;

    FileDimInfo out;
    ss >> out;

    EXPECT_EQ(out.name, "Intensity");
    EXPECT_EQ(out.type, pdal::Dimension::Type::Unsigned16);
    EXPECT_EQ(out.offset, 12);
}
