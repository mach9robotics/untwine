/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                      *
 *                                                                           *
 *   Unit tests for the chunker's re-chunk pass (chunker/ReChunk.cpp): the   *
 *   post-distribute recursive split that restores the per-chunk point       *
 *   budget when a count-grid cell over the budget became its own chunk.     *
 ****************************************************************************/

#include <gtest/gtest.h>

#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <vector>

#include "chunker/ChunkLut.hpp"
#include "chunker/Constants.hpp"
#include "chunker/ReChunk.hpp"
#include "untwine/Point.hpp"
#include "untwine/VoxelKey.hpp"

using namespace untwine;
using namespace untwine::chunker;

namespace
{

namespace fs = std::filesystem;

// A packed test point: x/y/z doubles (the on-disk layout Point reads) plus a payload the split
// must carry through byte-for-byte.
constexpr size_t PointSize = 3 * sizeof(double) + sizeof(uint64_t);

std::vector<char> makePoint(double x, double y, double z, uint64_t id)
{
    std::vector<char> rec(PointSize);
    std::memcpy(rec.data(), &x, sizeof(x));
    std::memcpy(rec.data() + 8, &y, sizeof(y));
    std::memcpy(rec.data() + 16, &z, sizeof(z));
    std::memcpy(rec.data() + 24, &id, sizeof(id));
    return rec;
}

struct TestPoint
{
    double x, y, z;
    uint64_t id;
};

void writeBin(const std::string& dir, const VoxelKey& key, const std::vector<TestPoint>& pts)
{
    std::ofstream out(dir + "/" + key.toString() + ".bin", std::ios::binary | std::ios::trunc);
    ASSERT_TRUE(out);
    for (const TestPoint& p : pts)
    {
        std::vector<char> rec = makePoint(p.x, p.y, p.z, p.id);
        out.write(rec.data(), rec.size());
    }
}

std::vector<TestPoint> readBin(const fs::path& path)
{
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in);
    std::vector<TestPoint> pts;
    std::vector<char> rec(PointSize);
    while (in.read(rec.data(), rec.size()))
    {
        TestPoint p;
        std::memcpy(&p.x, rec.data(), 8);
        std::memcpy(&p.y, rec.data() + 8, 8);
        std::memcpy(&p.z, rec.data() + 16, 8);
        std::memcpy(&p.id, rec.data() + 24, 8);
        pts.push_back(p);
    }
    return pts;
}

// All chunk .bins in `dir`, keyed by VoxelKey parsed from the filename.
std::map<std::string, std::vector<TestPoint>> readAllBins(const std::string& dir)
{
    std::map<std::string, std::vector<TestPoint>> bins;
    for (const auto& e : fs::directory_iterator(dir))
        if (e.path().extension() == ".bin")
            bins[e.path().stem().string()] = readBin(e.path());
    return bins;
}

// A fresh empty temp dir per test.
std::string testDir(const std::string& name)
{
    fs::path d = fs::path(::testing::TempDir()) / ("rechunk_" + name);
    fs::remove_all(d);
    fs::create_directories(d);
    return d.string();
}

pdal::BOX3D cube(double lo, double hi)
{
    return pdal::BOX3D(lo, lo, lo, hi, hi, hi);
}

} // namespace

// voxelBounds subdivides the full bounds by powers of two with the key selecting the cell —
// the same subdivision bu::VoxelInfo::bounds() computes.
TEST(ReChunk, voxel_bounds_subdivision)
{
    const pdal::BOX3D full = cube(0, 8);

    pdal::BOX3D b = voxelBounds(full, VoxelKey(0, 0, 0, 0));
    EXPECT_EQ(b, full);

    b = voxelBounds(full, VoxelKey(1, 0, 1, 1));
    EXPECT_DOUBLE_EQ(b.minx, 4); EXPECT_DOUBLE_EQ(b.maxx, 8);
    EXPECT_DOUBLE_EQ(b.miny, 0); EXPECT_DOUBLE_EQ(b.maxy, 4);
    EXPECT_DOUBLE_EQ(b.minz, 4); EXPECT_DOUBLE_EQ(b.maxz, 8);

    b = voxelBounds(full, VoxelKey(3, 0, 0, 2));
    EXPECT_DOUBLE_EQ(b.minx, 6); EXPECT_DOUBLE_EQ(b.maxx, 8);
    EXPECT_DOUBLE_EQ(b.miny, 0); EXPECT_DOUBLE_EQ(b.maxy, 2);
}

// countOctants bins by the node midpoint with the >= boundary rule and the VoxelKey::child bit
// layout (bit 0 x, bit 1 y, bit 2 z).
TEST(ReChunk, count_octants_midpoint_and_bit_layout)
{
    const pdal::BOX3D b = cube(0, 8);   // midpoint (4,4,4)

    std::vector<char> data;
    auto add = [&](double x, double y, double z)
    {
        std::vector<char> rec = makePoint(x, y, z, 0);
        data.insert(data.end(), rec.begin(), rec.end());
    };
    add(1, 1, 1);   // low corner -> octant 0
    add(5, 1, 1);   // +x -> octant 1
    add(1, 5, 1);   // +y -> octant 2
    add(1, 1, 5);   // +z -> octant 4
    add(5, 5, 5);   // high corner -> octant 7
    add(4, 4, 4);   // exactly on the midpoint -> upper side, octant 7

    std::array<uint64_t, 8> c = countOctants(data.data(), 6, PointSize, b);
    EXPECT_EQ(c[0], 1u);
    EXPECT_EQ(c[1], 1u);
    EXPECT_EQ(c[2], 1u);
    EXPECT_EQ(c[3], 0u);
    EXPECT_EQ(c[4], 1u);
    EXPECT_EQ(c[5], 0u);
    EXPECT_EQ(c[6], 0u);
    EXPECT_EQ(c[7], 2u);
}

// A chunk under the budget is left alone.
TEST(ReChunk, under_budget_chunk_untouched)
{
    const std::string dir = testDir("untouched");
    const VoxelKey key(0, 0, 0, 1);
    writeBin(dir, key, { {1, 1, 1, 10}, {2, 2, 2, 11} });

    std::vector<VoxelKey> out = rechunkFile(dir, cube(0, 8), PointSize, key, 5);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0], key);
    EXPECT_EQ(readAllBins(dir).size(), 1u);
    EXPECT_EQ(readAllBins(dir).at(key.toString()).size(), 2u);
}

// An oversized chunk splits into child .bins: the parent file is gone, every point survives with
// its payload intact, every point lands in the child whose bounds contain it, and every resulting
// chunk respects the budget.
TEST(ReChunk, oversized_chunk_splits_losslessly_under_budget)
{
    const std::string dir = testDir("split");
    const pdal::BOX3D full = cube(0, 8);
    const VoxelKey key(0, 0, 0, 1);   // bounds (0..4)^3, midpoint (2,2,2)

    // 12 points, 3 in each of 4 octants, ids 0..11.
    std::vector<TestPoint> pts;
    uint64_t id = 0;
    for (double d : {0.5, 1.0, 1.5})
    {
        pts.push_back({d, d, d, id++});            // octant 0
        pts.push_back({d + 2, d, d, id++});        // octant 1
        pts.push_back({d, d + 2, d, id++});        // octant 2
        pts.push_back({d + 2, d + 2, d + 2, id++}); // octant 7
    }
    writeBin(dir, key, pts);

    std::vector<VoxelKey> out = rechunkFile(dir, full, PointSize, key, 4);

    auto bins = readAllBins(dir);
    EXPECT_EQ(bins.count(key.toString()), 0u) << "parent .bin should be deleted";
    EXPECT_EQ(bins.size(), out.size());

    size_t total = 0;
    std::set<uint64_t> ids;
    for (const VoxelKey& k : out)
    {
        auto it = bins.find(k.toString());
        ASSERT_NE(it, bins.end()) << "missing .bin for returned chunk " << k;
        EXPECT_LE(it->second.size(), 4u) << "chunk " << k << " over budget";
        const pdal::BOX3D b = voxelBounds(full, k);
        for (const TestPoint& p : it->second)
        {
            EXPECT_TRUE(b.contains(p.x, p.y, p.z))
                << "point " << p.id << " outside chunk " << k;
            ids.insert(p.id);
            ++total;
        }
    }
    EXPECT_EQ(total, pts.size());
    EXPECT_EQ(ids.size(), pts.size()) << "ids must survive the split uniquely";
}

// A tight cluster occupying one octant descends by rename (no rewrite) until it separates, so the
// resulting chunks sit several levels below the original key.
TEST(ReChunk, single_octant_cluster_descends_until_it_separates)
{
    const std::string dir = testDir("descend");
    const pdal::BOX3D full = cube(0, 8);
    const VoxelKey key(0, 0, 0, 1);   // bounds (0..4)^3

    // 6 points inside (0..0.5)^3: one occupied octant at levels 1, 2, and 3. They separate at
    // (0..0.5) midpoint 0.25.
    std::vector<TestPoint> pts;
    for (int i = 0; i < 3; ++i)
        pts.push_back({0.1, 0.1, 0.1 + i * 0.01, (uint64_t)i});
    for (int i = 0; i < 3; ++i)
        pts.push_back({0.4, 0.4, 0.4 + i * 0.01, (uint64_t)(10 + i)});
    writeBin(dir, key, pts);

    std::vector<VoxelKey> out = rechunkFile(dir, full, PointSize, key, 4);

    ASSERT_EQ(out.size(), 2u);
    for (const VoxelKey& k : out)
        EXPECT_GE(k.level(), 4) << "cluster should have descended below level 3, got " << k;

    auto bins = readAllBins(dir);
    EXPECT_EQ(bins.size(), 2u);
    size_t total = 0;
    for (const auto& b : bins)
    {
        EXPECT_LE(b.second.size(), 4u);
        total += b.second.size();
    }
    EXPECT_EQ(total, pts.size());
}

// Coincident points can never separate: recursion must stop at MaxRechunkLevel and leave one
// (still oversized) chunk rather than loop forever. The Indexing phase guards this residual.
TEST(ReChunk, coincident_points_stop_at_level_cap)
{
    const std::string dir = testDir("coincident");
    const VoxelKey key(0, 0, 0, 1);

    std::vector<TestPoint> pts(10, TestPoint{1.5, 1.5, 1.5, 0});
    for (size_t i = 0; i < pts.size(); ++i)
        pts[i].id = i;
    writeBin(dir, key, pts);

    std::vector<VoxelKey> out = rechunkFile(dir, cube(0, 8), PointSize, key, 4);

    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].level(), MaxRechunkLevel);
    auto bins = readAllBins(dir);
    ASSERT_EQ(bins.size(), 1u);
    EXPECT_EQ(bins.at(out[0].toString()).size(), pts.size());
}

// rechunkOversized touches only chunks whose planned count exceeds the target: the oversized
// cell's chunk is split, the under-budget chunk's .bin is left byte-identical.
TEST(ReChunk, only_oversized_chunks_are_split)
{
    const std::string dir = testDir("oversized_only");
    const pdal::BOX3D full = cube(0, 8);

    const VoxelKey small(0, 0, 0, 2);   // bounds (0..2)^3
    const VoxelKey big(1, 1, 1, 2);     // bounds (2..4)^3

    std::vector<TestPoint> smallPts{ {0.5, 0.5, 0.5, 100}, {1.5, 1.5, 1.5, 101} };
    std::vector<TestPoint> bigPts;
    for (int i = 0; i < 8; ++i)   // one point per octant of (2..4)^3
        bigPts.push_back({2.5 + (i & 1), 2.5 + ((i >> 1) & 1), 2.5 + ((i >> 2) & 1),
            (uint64_t)i});
    writeBin(dir, small, smallPts);
    writeBin(dir, big, bigPts);

    ChunkLut lut;
    lut.roots = { small, big };
    lut.cellToRoot = { { small, small }, { big, big } };
    std::unordered_map<VoxelKey, uint64_t> cellCounts{
        { small, smallPts.size() }, { big, bigPts.size() } };

    const size_t split = rechunkOversized(dir, full, PointSize, lut, cellCounts, 4);
    EXPECT_EQ(split, 1u);

    auto bins = readAllBins(dir);
    ASSERT_TRUE(bins.count(small.toString())) << "under-budget chunk must be untouched";
    EXPECT_EQ(bins.at(small.toString()).size(), smallPts.size());
    EXPECT_EQ(bins.count(big.toString()), 0u) << "oversized chunk should have been split";

    size_t total = 0;
    for (const auto& b : bins)
    {
        EXPECT_LE(b.second.size(), 4u);
        total += b.second.size();
    }
    EXPECT_EQ(total, smallPts.size() + bigPts.size());
}
