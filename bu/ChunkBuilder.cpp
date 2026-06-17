/*****************************************************************************
 *   Copyright (c) 2020, Hobu, Inc. (info@hobu.co)                           *
 *   Modified by Kyle Kam (kylekam@mach9.io)                                 *
 *                                                                           *
 *   All rights reserved.                                                    *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 ****************************************************************************/

#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <numeric>
#include <thread>
#include <unordered_set>
#include <vector>

#include <pdal/PDALUtils.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/filters/SortFilter.hpp>
#include <pdal/util/Algorithm.hpp>
#include <pdal/util/FileUtils.hpp>

#include <lazperf/lazperf.hpp>
#include <lazperf/writers.hpp>
#include <lazperf/readers.hpp>

#include "../untwine/Common.hpp"
#include "../untwine/GridKey.hpp"
#include "../untwine/Point.hpp"
#include "../untwine/ProgressWriter.hpp"
#include "../untwine/ThreadPool.hpp"
#include "../untwine/VoxelKey.hpp"

#include "ChunkBuilder.hpp"
#include "PointAccessor.hpp"
#include "PyramidManager.hpp"
#include "Stats.hpp"
#include "VoxelInfo.hpp"

#include <stringconv.hpp>  // untwine/os

namespace untwine
{
namespace bu
{

namespace
{

// One chunk to build: an internal subtree root plus its descendant leaf files.
struct ChunkPlan
{
    VoxelKey root;
    std::unordered_map<VoxelKey, FileInfo> leaves;
};

// Partition the leaves into adaptively-sized internal chunk roots.
//
// Build a per-key subtree point count, then descend from the global root: stop at the
// shallowest key whose subtree holds <= target points. A stopping key that is itself a leaf
// is a "single-leaf chunk root" -- left on disk untouched for the cap, not returned here.
// Internal stopping keys become ChunkPlans. The global root never stops (chunk roots are at
// level >= 1), so the cap always has at least the root to process.
std::vector<ChunkPlan> planChunks(const std::unordered_map<VoxelKey, FileInfo>& leaves,
    uint64_t target)
{
    const VoxelKey root;

    std::unordered_map<VoxelKey, uint64_t> cnt;
    for (const auto& p : leaves)
    {
        VoxelKey k = p.first;
        const uint64_t n = (uint64_t)p.second.numPoints();
        while (true)
        {
            cnt[k] += n;
            if (k == root)
                break;
            k = k.parent();
        }
    }

    std::vector<VoxelKey> internalRoots;
    std::function<void(const VoxelKey&)> descend = [&](const VoxelKey& k)
    {
        auto ci = cnt.find(k);
        if (ci == cnt.end() || ci->second == 0)
            return;
        bool isLeaf = leaves.count(k) > 0;
        bool fits = ci->second <= target;
        if (k != root && (isLeaf || fits))
        {
            if (!isLeaf)
                internalRoots.push_back(k);
            return;
        }
        for (int i = 0; i < 8; ++i)
            descend(k.child(i));
    };
    descend(root);

    std::unordered_map<VoxelKey, ChunkPlan> byRoot;
    for (const VoxelKey& r : internalRoots)
        byRoot[r].root = r;

    // Assign each leaf to its nearest internal-root ancestor (if any).
    for (const auto& p : leaves)
    {
        VoxelKey k = p.first;
        while (true)
        {
            auto bi = byRoot.find(k);
            if (bi != byRoot.end())
            {
                bi->second.leaves.emplace(p.first, p.second);
                break;
            }
            if (k == root)
                break;
            k = k.parent();
        }
    }

    std::vector<ChunkPlan> plans;
    plans.reserve(byRoot.size());
    for (auto& bp : byRoot)
        plans.push_back(std::move(bp.second));
    return plans;
}

// Builds one chunk subtree fully in RAM and emits COPC chunks for its descendant nodes.
class ChunkBuilder
{
public:
    ChunkBuilder(const BaseInfo& b, PyramidManager& mgr, ChunkPlan&& plan) :
        m_b(b), m_mgr(mgr), m_root(plan.root), m_leaves(std::move(plan.leaves)), m_points(b)
    {
        for (const auto& p : m_leaves)
        {
            VoxelKey k = p.first;
            while (true)
            {
                m_present.insert(k);
                if (k == m_root)
                    break;
                k = k.parent();
            }
        }
    }

    void run()
    {
        std::vector<int> accepted = build(m_root);
        writeBin(m_root, accepted);
    }

private:
    // Returns the indices (into m_points) accepted at (promoted to) node k, emitting a COPC
    // chunk for each of k's children along the way.
    std::vector<int> build(const VoxelKey& k);
    void emit(const VoxelKey& key, const std::vector<int>& indices);
    void createChunk(const VoxelKey& key, pdal::PointViewPtr view,
        const pdal::DimTypeList& extraDims);
    void fillPointBuf(pdal::PointRef& point, std::vector<char>& buf,
        pdal::Dimension::Id bitsDim, const pdal::DimTypeList& extraDims);
    void sortChunk(pdal::PointViewPtr view);
    void writeBin(const VoxelKey& key, const std::vector<int>& indices);

    const BaseInfo& m_b;
    PyramidManager& m_mgr;
    VoxelKey m_root;
    std::unordered_map<VoxelKey, FileInfo> m_leaves;
    std::unordered_set<VoxelKey> m_present;
    PointAccessor m_points;
};

std::vector<int> ChunkBuilder::build(const VoxelKey& k)
{
    // Leaf: map its .bin into the accessor and hand up its full index range.
    auto li = m_leaves.find(k);
    if (li != m_leaves.end())
    {
        FileInfo& fi = li->second;
        int start = (int)m_points.size();
        m_points.read(fi);
        std::vector<int> idx(fi.numPoints());
        std::iota(idx.begin(), idx.end(), start);
        return idx;
    }

    std::array<std::vector<int>, 8> childAccepted;
    for (int dir = 0; dir < 8; ++dir)
    {
        VoxelKey c = k.child(dir);
        if (m_present.count(c))
            childAccepted[dir] = build(c);
    }

    // Sample one point per occupied cell of k's grid (sequential order, matching the
    // de-randomized default sampler). Accepted points are promoted to k; the rest are
    // rejected back to the child they came from and become that child's COPC chunk.
    VoxelInfo vi(m_b.bounds, k);
    VoxelInfo::Grid& grid = vi.grid();
    std::vector<int> accepted;
    std::array<std::vector<int>, 8> rejected;
    for (int dir = 0; dir < 8; ++dir)
    {
        for (int idx : childAccepted[dir])
        {
            const Point& p = m_points[idx];
            GridKey gk = vi.gridKey(p);
            if (grid.find(gk) == grid.end())
            {
                grid.insert({gk, idx});
                accepted.push_back(idx);
            }
            else
                rejected[dir].push_back(idx);
        }
    }

    // Any child that contributed points must appear in the tree; emit its chunk (the
    // rejected list may be empty -> a structural node with zero points but children).
    for (int dir = 0; dir < 8; ++dir)
        if (!childAccepted[dir].empty())
            emit(k.child(dir), rejected[dir]);

    return accepted;
}

void ChunkBuilder::emit(const VoxelKey& key, const std::vector<int>& indices)
{
    using namespace pdal;

    PointTable table;
    IndexedStats stats;
    DimTypeList extraDims;
    DimInfoList dims = m_b.dimInfo;
    for (FileDimInfo& fdi : dims)
    {
        fdi.dim = table.layout()->registerOrAssignDim(fdi.name, fdi.type);
        if (m_b.opts.stats)
        {
            if (fdi.dim == Dimension::Id::Classification)
                stats.push_back({fdi.dim, Stats(fdi.name, Stats::EnumType::Enumerate, false)});
            else if (fdi.dim == Dimension::Id::ReturnNumber)
                stats.push_back({fdi.dim, Stats(fdi.name, Stats::EnumType::Enumerate, false)});
            else
                stats.push_back({fdi.dim, Stats(fdi.name, Stats::EnumType::NoEnum, false)});
        }
        if (fdi.extraDim)
            extraDims.push_back(DimType(fdi.dim, fdi.type));
    }
    table.finalize();

    PointViewPtr view(new PointView(table));
    PointId pointId = 0;
    for (int idx : indices)
    {
        char *base = m_points[idx].cdata();
        for (const FileDimInfo& fdi : dims)
            view->setField(fdi.dim, fdi.type, pointId,
                reinterpret_cast<void *>(base + fdi.offset));
        pointId++;
    }

    if (m_b.opts.stats)
    {
        for (PointId id = 0; id < view->size(); ++id)
            for (auto& sp : stats)
                sp.second.insert(view->getFieldAs<double>(sp.first, id));
    }

    try
    {
        createChunk(key, view, extraDims);
    }
    catch (pdal_error& err)
    {
        throw FatalError(err.what());
    }
    m_mgr.logOctant(key, (int)indices.size(), stats);
}

void ChunkBuilder::createChunk(const VoxelKey& key, pdal::PointViewPtr view,
    const pdal::DimTypeList& extraDims)
{
    using namespace pdal;

    if (view->size() == 0)
    {
        m_mgr.newChunk(key, 0, 0);
        return;
    }

    if (view->layout()->hasDim(Dimension::Id::GpsTime))
        sortChunk(view);

    PointLayoutPtr layout = view->layout();

    int ebCount {0};
    for (DimType dim : extraDims)
        ebCount += layout->dimSize(dim.m_id);

    std::vector<char> buf(lazperf::baseCount(m_b.pointFormatId) + ebCount);
    lazperf::writer::chunk_compressor compressor(m_b.pointFormatId, ebCount);
    for (PointId idx = 0; idx < view->size(); ++idx)
    {
        PointRef point(*view, idx);
        fillPointBuf(point, buf, layout->findDim(UntwineBitsDimName), extraDims);
        compressor.compress(buf.data());
    }
    std::vector<unsigned char> chunk = compressor.done();

    uint64_t location = m_mgr.newChunk(key, chunk.size(), (uint32_t)view->size());

    std::ofstream out(os::toNative(m_b.opts.outputName),
        std::ios::out | std::ios::in | std::ios::binary);
    out.seekp(std::ofstream::pos_type(location));
    out.write(reinterpret_cast<const char *>(chunk.data()), chunk.size());
    out.close();
    if (!out)
        throw FatalError("Failure writing to file '" + m_b.opts.outputName + "'.");
}

void ChunkBuilder::sortChunk(pdal::PointViewPtr view)
{
    pdal::BufferReader r;
    r.addView(view);

    pdal::SortFilter s;
    s.setInput(r);
    pdal::Options o;
    o.add("dimension", "GpsTime");
    s.setOptions(o);

    s.prepare(view->table());
    s.execute(view->table());
}

void ChunkBuilder::fillPointBuf(pdal::PointRef& point, std::vector<char>& buf,
    pdal::Dimension::Id bitsDim, const pdal::DimTypeList& extraDims)
{
    using namespace pdal;

    LeInserter ostream(buf.data(), buf.size());

    bool hasTime = true;
    bool hasColor = m_b.pointFormatId == 7 || m_b.pointFormatId == 8;
    bool hasInfrared = m_b.pointFormatId == 8;

    using namespace Dimension;

    uint8_t returnNumber(1);
    uint8_t numberOfReturns(1);
    if (point.hasDim(Id::ReturnNumber))
        returnNumber = point.getFieldAs<uint8_t>(Id::ReturnNumber);
    if (point.hasDim(Id::NumberOfReturns))
        numberOfReturns = point.getFieldAs<uint8_t>(Id::NumberOfReturns);

    auto converter = [](double d, Dimension::Id dim) -> int32_t
    {
        int32_t i(0);

        if (!Utils::numericCast(d, i))
            throw FatalError("Unable to convert scaled value (" +
                Utils::toString(d) + ") to "
                "int32 for dimension '" + Dimension::name(dim) +
                "' when writing LAS/LAZ file.");
        return i;
    };

    double x = point.getFieldAs<double>(Id::X);
    int32_t xi = converter((x - m_b.xform.offset.x) / m_b.xform.scale.x, Id::X);
    double y = point.getFieldAs<double>(Id::Y);
    int32_t yi = converter((y - m_b.xform.offset.y) / m_b.xform.scale.y, Id::Y);
    double z = point.getFieldAs<double>(Id::Z);
    int32_t zi = converter((z - m_b.xform.offset.z) / m_b.xform.scale.z, Id::Z);

    ostream << xi << yi << zi;

    ostream << point.getFieldAs<uint16_t>(Id::Intensity);
    ostream << (uint8_t)(returnNumber | (numberOfReturns << 4));
    ostream << point.getFieldAs<uint8_t>(bitsDim);
    ostream << point.getFieldAs<uint8_t>(Id::Classification);

    uint8_t userData = point.getFieldAs<uint8_t>(Id::UserData);
    int16_t scanAngleRank =
        static_cast<int16_t>(std::round(
            point.getFieldAs<float>(Id::ScanAngleRank) / .006f));
    ostream << userData << scanAngleRank;

    ostream << point.getFieldAs<uint16_t>(Id::PointSourceId);

    if (hasTime)
        ostream << point.getFieldAs<double>(Id::GpsTime);

    if (hasColor)
    {
        ostream << point.getFieldAs<uint16_t>(Id::Red);
        ostream << point.getFieldAs<uint16_t>(Id::Green);
        ostream << point.getFieldAs<uint16_t>(Id::Blue);
    }

    if (hasInfrared)
        ostream << point.getFieldAs<uint16_t>(Id::Infrared);

    Everything e;
    for (auto& dim : extraDims)
    {
        point.getField((char *)&e, dim.m_id, dim.m_type);
        Utils::insertDim(ostream, dim.m_type, e);
    }
}

void ChunkBuilder::writeBin(const VoxelKey& key, const std::vector<int>& indices)
{
    if (indices.empty())
        return;

    std::string fullFilename = m_b.opts.tempDir + "/" + key.toString() + ".bin";
    std::ofstream out(os::toNative(fullFilename), std::ios::binary | std::ios::trunc);
    if (!out)
        throw FatalError("Couldn't open '" + fullFilename + "' for output.");
    for (int i : indices)
        out.write(m_points[i].cdata(), m_b.pointSize);
}

} // unnamed namespace

void runChunkedBuild(const BaseInfo& b, PyramidManager& mgr,
    const std::unordered_map<VoxelKey, FileInfo>& leaves,
    uint64_t maxChunkPoints, ProgressWriter* progress)
{
    std::vector<ChunkPlan> plans = planChunks(leaves, maxChunkPoints);
    {
        uint64_t chunkLeaves = 0;
        for (const ChunkPlan& p : plans)
            chunkLeaves += p.leaves.size();
        std::cerr << "[chunked-build] target=" << maxChunkPoints << " leaves=" << leaves.size()
                  << " internal-chunk-roots=" << plans.size()
                  << " leaves-consumed=" << chunkLeaves << "\n";
    }
    if (plans.empty())
        return;

    // Per-chunk increment: advance the bar 60->90% over the build (cap finishes 90->100%).
    if (progress)
        progress->setIncrement(0.30 / (double)plans.size());

    size_t nthreads = std::thread::hardware_concurrency();
    if (nthreads == 0)
        nthreads = 8;

    ThreadPool pool(nthreads);
    pool.trap(true);

    for (ChunkPlan& plan : plans)
    {
        ChunkPlan *pp = &plan;
        pool.add([pp, &b, &mgr, progress]()
        {
            // Filenames to delete once the builder has consumed (and unmapped) them.
            std::vector<std::string> consumed;
            consumed.reserve(pp->leaves.size());
            for (const auto& l : pp->leaves)
                consumed.push_back(l.second.filename());

            {
                ChunkBuilder cb(b, mgr, std::move(*pp));
                cb.run();
            } // ChunkBuilder/PointAccessor destroyed here -> leaf files unmapped

            for (const std::string& fn : consumed)
                pdal::FileUtils::deleteFile(b.opts.tempDir + "/" + fn);
            if (progress)
                progress->writeIncrement("Chunk built");
        });
    }

    pool.await();
    if (pool.hasErrors())
    {
        std::vector<std::string> errs = pool.clearErrors();
        throw FatalError(errs.empty() ? std::string("Chunked build failed.") : errs.front());
    }
    pool.join();
}

} // namespace bu
} // namespace untwine
