/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                       *
 *                                                                           *
 *   Unit tests for bu::Stats -- the per-dimension summary the BU step        *
 *   accumulates per node and merges up the tree into the whole-cloud stats   *
 *   written to the COPC header. merge() correctness is a BU invariant.      *
 ****************************************************************************/

#include <gtest/gtest.h>

#include "bu/Stats.hpp"

using namespace untwine;

namespace
{
const std::vector<double> kSeq { 2, 4, 4, 4, 5, 5, 7, 9 };  // mean 5, pop var 4
}

// A single insert: min == max == average == the value, count 1.
TEST(Stats, single_value)
{
    Stats s("X", Stats::NoEnum);
    s.insert(42.0);
    EXPECT_EQ(s.count(), 1u);
    EXPECT_DOUBLE_EQ(s.minimum(), 42.0);
    EXPECT_DOUBLE_EQ(s.maximum(), 42.0);
    EXPECT_DOUBLE_EQ(s.average(), 42.0);
}

// Known sequence -> known moments (textbook population stddev of 2).
TEST(Stats, known_moments)
{
    Stats s("X", Stats::NoEnum);
    for (double v : kSeq)
        s.insert(v);

    EXPECT_EQ(s.count(), kSeq.size());
    EXPECT_DOUBLE_EQ(s.minimum(), 2.0);
    EXPECT_DOUBLE_EQ(s.maximum(), 9.0);
    EXPECT_NEAR(s.average(), 5.0, 1e-12);
    EXPECT_NEAR(s.populationVariance(), 4.0, 1e-12);
    EXPECT_NEAR(s.populationStddev(), 2.0, 1e-12);
}

// THE BU rollup invariant: accumulating per node then merging must equal
// accumulating the whole cloud in one pass -- count, extrema and moments.
TEST(Stats, merge_equals_single_pass)
{
    Stats whole("X", Stats::NoEnum);
    for (double v : kSeq)
        whole.insert(v);

    // Split the sequence across two "nodes" and merge.
    Stats a("X", Stats::NoEnum);
    Stats b("X", Stats::NoEnum);
    for (size_t i = 0; i < kSeq.size(); ++i)
        (i < 3 ? a : b).insert(kSeq[i]);

    ASSERT_TRUE(a.merge(b));

    EXPECT_EQ(a.count(), whole.count());
    EXPECT_DOUBLE_EQ(a.minimum(), whole.minimum());
    EXPECT_DOUBLE_EQ(a.maximum(), whole.maximum());
    EXPECT_NEAR(a.average(), whole.average(), 1e-9);
    EXPECT_NEAR(a.populationVariance(), whole.populationVariance(), 1e-9);
    EXPECT_NEAR(a.skewness(), whole.skewness(), 1e-9);
    EXPECT_NEAR(a.kurtosis(), whole.kurtosis(), 1e-9);
}

// merge refuses mismatched summaries (different name / enum / advanced flag),
// guarding against accidentally folding the wrong dimension together.
TEST(Stats, merge_rejects_mismatch)
{
    Stats x("X", Stats::NoEnum);
    Stats y("Y", Stats::NoEnum);
    x.insert(1.0);
    y.insert(2.0);
    EXPECT_FALSE(x.merge(y));

    Stats enumX("X", Stats::Enumerate);
    EXPECT_FALSE(x.merge(enumX));
}

// Enumeration mode tracks value -> occurrence counts.
TEST(Stats, enumerate_counts)
{
    Stats s("Classification", Stats::Enumerate);
    for (double v : { 1.0, 1.0, 2.0, 1.0, 2.0 })
        s.insert(v);

    const auto& m = s.values();
    EXPECT_EQ(m.at(1.0), 3u);
    EXPECT_EQ(m.at(2.0), 2u);
}

// Global mode retains data so median + MAD can be computed.
TEST(Stats, global_median_and_mad)
{
    Stats s("X", Stats::Global);
    for (double v : { 1.0, 2.0, 3.0, 4.0, 5.0 })
        s.insert(v);
    s.computeGlobalStats();

    EXPECT_DOUBLE_EQ(s.median(), 3.0);
    // abs deviations from median {2,1,0,1,2} -> median 1.
    EXPECT_DOUBLE_EQ(s.mad(), 1.0);
}
