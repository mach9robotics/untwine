/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Unit tests for Point (the xyz view over a packed record) and Transform   *
 *   (scale/offset validity). Point::quantize is applied during EPF so that   *
 *   binned coordinates snap to the output scale grid.                       *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <array>
#include <cstring>

#include "untwine/Point.hpp"
#include "untwine/Common.hpp"

using namespace untwine;

namespace
{
// Build a packed record whose first three doubles are x, y, z.
std::array<char, 64> makeRecord(double x, double y, double z)
{
    std::array<char, 64> buf{};
    std::memcpy(buf.data() + 0 * sizeof(double), &x, sizeof(double));
    std::memcpy(buf.data() + 1 * sizeof(double), &y, sizeof(double));
    std::memcpy(buf.data() + 2 * sizeof(double), &z, sizeof(double));
    return buf;
}
}

// x()/y()/z() read the leading three doubles of the record.
TEST(Point, reads_xyz)
{
    auto rec = makeRecord(1.5, -2.25, 100.125);
    Point p(rec.data());
    EXPECT_DOUBLE_EQ(p.x(), 1.5);
    EXPECT_DOUBLE_EQ(p.y(), -2.25);
    EXPECT_DOUBLE_EQ(p.z(), 100.125);
}

// quantize snaps each coordinate to round((v-offset)/scale)*scale+offset.
TEST(Point, quantize_snaps_to_scale_grid)
{
    auto rec = makeRecord(1.2345, 5.6789, -9.8765);
    Point p(rec.data());

    Transform xform;
    xform.scale = { 0.01, 0.01, 0.01 };
    xform.offset = { 0.0, 0.0, 0.0 };
    p.quantize(xform);

    EXPECT_NEAR(p.x(), 1.23, 1e-9);
    EXPECT_NEAR(p.y(), 5.68, 1e-9);
    EXPECT_NEAR(p.z(), -9.88, 1e-9);
}

// quantize respects a non-zero offset (the grid is anchored at offset).
TEST(Point, quantize_with_offset)
{
    auto rec = makeRecord(10.04, 0, 0);
    Point p(rec.data());

    Transform xform;
    xform.scale = { 0.1, 0.1, 0.1 };
    xform.offset = { 10.0, 0.0, 0.0 };
    p.quantize(xform);

    // (10.04 - 10.0)/0.1 = 0.4 -> round 0 -> 10.0
    EXPECT_NEAR(p.x(), 10.0, 1e-9);
}

// A default Transform is invalid (zero scale, NaN offset); a fully specified
// one is valid. This gates whether prep computed a usable transform.
TEST(Transform, validity)
{
    Transform def;
    EXPECT_FALSE(def.valid());

    Transform zeroScale;
    zeroScale.scale = { 0.0, 1.0, 1.0 };
    zeroScale.offset = { 0.0, 0.0, 0.0 };
    EXPECT_FALSE(zeroScale.valid());

    Transform good;
    good.scale = { 0.01, 0.01, 0.01 };
    good.offset = { 0.0, 0.0, 0.0 };
    EXPECT_TRUE(good.valid());
}
