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
#include "ChunkPlan.hpp"
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

// Leaf cutoff for the in-chunk octree, matching epf MaxPointsPerNode, plus a depth guard so
// pathological coincident points can't recurse forever.
constexpr size_t LeafPointBudget = 100000;
constexpr int MaxBuildLevel = 30;

// Fallback worker-thread count when std::thread::hardware_concurrency() reports 0.
constexpr unsigned DefaultThreadCount = 8;

// One chunk to build: an internal subtree root plus its descendant leaf files.
struct ChunkPlan
{
    VoxelKey root;
    std::unordered_map<VoxelKey, FileInfo> leaves;
};

// Partition the leaves into adaptively-sized internal chunk roots.
//
// Build a per-key subtree point count, then descend from the global root and stop at the
// shallowest key whose subtree holds <= target points. A stopping key that is itself a leaf is a
// single-leaf chunk root, left on disk untouched for the merge and not returned here. Internal
// stopping keys become ChunkPlans. The global root never stops, so chunk roots are at level >= 1
// and the merge always has at least the root to process.
std::vector<ChunkPlan> planChunks(const std::unordered_map<VoxelKey, FileInfo>& leaves,
    uint64_t target)
{
    const VoxelKey root;

    // Internal chunk roots from the leaf point counts, via the pure logic in ChunkPlan.
    std::unordered_map<VoxelKey, uint64_t> counts;
    counts.reserve(leaves.size());
    for (const auto& p : leaves)
        counts[p.first] = (uint64_t)p.second.numPoints();

    std::unordered_map<VoxelKey, ChunkPlan> byRoot;
    for (const VoxelKey& r : planChunkRoots(counts, target))
        byRoot[r].root = r;

    // Assign each leaf to its nearest internal-root ancestor, if any.
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

// Indexing: builds one chunk's local octree in RAM and emits a COPC chunk per node.
class ChunkBuilder
{
public:
    ChunkBuilder(const BaseInfo& b, PyramidManager& mgr, ChunkPlan&& plan) :
        m_b(b), m_mgr(mgr), m_root(plan.root), m_leaves(std::move(plan.leaves)), m_points(b)
    {}

    void run()
    {
        // Load the whole chunk's points into one flat index space.
        for (auto& p : m_leaves)
            m_points.read(p.second);
        std::vector<int> all(m_points.size());
        std::iota(all.begin(), all.end(), 0);

        // Phase 1: decide the chunk's octree structure by histogramming points into a per-chunk
        // local count grid and merging cells under the leaf budget, recursing into any cell still
        // over budget. This is PotreeConverter's per-chunk grid plus merge-under-count.
        buildStructure(m_root, std::move(all));

        // Phase 2: sample the structure bottom-up, one point per 128^3 cell, emitting each
        // descendant node's COPC chunk. The chunk root's promoted points go to <root>.bin.
        std::vector<int> accepted = sampleStructure(m_root);
        writeBin(m_root, accepted);
    }

private:
    // Phase 1. Recursively histogram `indices` into a local grid below `node`, merge
    // cells under LeafPointBudget into this node's octree leaves, and recurse into oversized cells.
    // Populates m_leafSet / m_leafPoints / m_occupied.
    void buildStructure(const VoxelKey& node, std::vector<int>&& indices);
    void markLeaf(const VoxelKey& node, std::vector<int>&& indices);
    // The depth-`depth` descendant of `startKey`, whose bounds are `startBounds`, holding point p.
    VoxelKey cellOf(const Point& p, const VoxelKey& startKey, const pdal::BOX3D& startBounds,
        int depth);
    // Phase 2. Sample node k bottom-up over the precomputed structure, emit children's chunks,
    // and return the points promoted to k.
    std::vector<int> sampleStructure(const VoxelKey& k);

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
    PointAccessor m_points;
    // Octree structure for this chunk, filled by buildStructure and consumed by sampleStructure.
    std::unordered_set<VoxelKey> m_leafSet;                       // the leaves
    std::unordered_map<VoxelKey, std::vector<int>> m_leafPoints;  // leaf -> its point indices
    std::unordered_set<VoxelKey> m_occupied;                      // every node that exists
};

VoxelKey ChunkBuilder::cellOf(const Point& p, const VoxelKey& startKey,
    const pdal::BOX3D& startBounds, int depth)
{
    VoxelKey k = startKey;
    pdal::BOX3D b = startBounds;
    for (int d = 0; d < depth; ++d)
    {
        const int oct = childOctant(p.x(), p.y(), p.z(), b);
        const double mx = (b.minx + b.maxx) * 0.5;
        const double my = (b.miny + b.maxy) * 0.5;
        const double mz = (b.minz + b.maxz) * 0.5;
        if (oct & 1) b.minx = mx; else b.maxx = mx;   // bit0 -> x
        if (oct & 2) b.miny = my; else b.maxy = my;   // bit1 -> y
        if (oct & 4) b.minz = mz; else b.maxz = mz;   // bit2 -> z
        k = k.child(oct);
    }
    return k;
}

void ChunkBuilder::markLeaf(const VoxelKey& node, std::vector<int>&& indices)
{
    m_leafSet.insert(node);
    m_leafPoints[node] = std::move(indices);
    // Mark the leaf and its ancestors up to the chunk root as existing nodes. Stop early once we
    // reach a node already marked, since its ancestors are then marked too.
    for (VoxelKey k = node;; k = k.parent())
    {
        if (!m_occupied.insert(k).second)
            break;
        if (k == m_root)
            break;
    }
}

void ChunkBuilder::buildStructure(const VoxelKey& node, std::vector<int>&& indices)
{
    // Small enough or too deep to split: this node is a leaf holding all its points.
    if (indices.size() <= LeafPointBudget || node.level() >= MaxBuildLevel)
    {
        markLeaf(node, std::move(indices));
        return;
    }

    // Histogram the points into a per-chunk local grid `depth` levels below `node`.
    const int depth = localGridDepth(indices.size(), LeafPointBudget,
        MaxBuildLevel - node.level());
    const pdal::BOX3D nodeBounds = VoxelInfo(m_b.bounds, node).bounds();

    std::unordered_map<VoxelKey, std::vector<int>> cells;
    for (int idx : indices)
        cells[cellOf(m_points[idx], node, nodeBounds, depth)].push_back(idx);

    // Merge cells under the leaf budget into the cut nodes, this node's octree leaves at this
    // scale. Same operation as the global chunk planner. A cell over budget comes back as its own
    // oversized cut node, which we then refine with a deeper grid.
    std::unordered_map<VoxelKey, uint64_t> counts;
    counts.reserve(cells.size());
    for (const auto& c : cells)
        counts[c.first] = (uint64_t)c.second.size();
    std::vector<VoxelKey> cut = planChunkRoots(counts, LeafPointBudget);
    std::unordered_set<VoxelKey> cutSet(cut.begin(), cut.end());

    // Group each cell's points under its covering cut node. As a fallback, a small standalone cell
    // that planChunkRoots left out becomes its own leaf.
    std::unordered_map<VoxelKey, std::vector<int>> cutPoints;
    for (auto& c : cells)
    {
        VoxelKey r = c.first;
        while (!cutSet.count(r) && r != node)
            r = r.parent();
        if (!cutSet.count(r))
            r = c.first;
        std::vector<int>& dst = cutPoints[r];
        dst.insert(dst.end(), c.second.begin(), c.second.end());
    }

    // A cut node still over budget is an oversized finest cell, so recurse with a fresh deeper
    // grid. Everything else is a finished leaf.
    for (auto& cp : cutPoints)
    {
        if (cp.second.size() > LeafPointBudget && cp.first.level() < MaxBuildLevel)
            buildStructure(cp.first, std::move(cp.second));
        else
            markLeaf(cp.first, std::move(cp.second));
    }
}

std::vector<int> ChunkBuilder::sampleStructure(const VoxelKey& k)
{
    // Leaf: hand its points up unchanged; the parent samples them and emits this node's chunk.
    if (m_leafSet.count(k))
        return std::move(m_leafPoints[k]);

    // Internal node: gather the points each occupied child promoted upward.
    std::array<std::vector<int>, 8> childAccepted;
    for (int dir = 0; dir < 8; ++dir)
    {
        const VoxelKey c = k.child(dir);
        if (m_occupied.count(c))
            childAccepted[dir] = sampleStructure(c);
    }

    // Sample one point per occupied cell of k's 128^3 grid. Accepted points are promoted to k;
    // the rest are rejected back to their child and become that child's COPC chunk. Order is
    // irrelevant, since only cell occupancy is tested.
    VoxelInfo vi(m_b.bounds, k);
    VoxelInfo::Grid& grid = vi.grid();
    std::vector<int> accepted;
    std::array<std::vector<int>, 8> rejected;
    for (int dir = 0; dir < 8; ++dir)
    {
        for (int idx : childAccepted[dir])
        {
            GridKey gk = vi.gridKey(m_points[idx]);
            if (grid.find(gk) == grid.end())
            {
                grid.insert({gk, idx});
                accepted.push_back(idx);
            }
            else
                rejected[dir].push_back(idx);
        }
    }

    // Any child that received points must appear in the tree, so emit its chunk. The rejected list
    // may be empty, giving a structural node with children but zero points.
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

    const std::string fullFilename = m_b.opts.tempDir + "/" + key.toString() + ".bin";
    const std::string tmpFilename = fullFilename + ".tmp";
    {
        std::ofstream out(os::toNative(tmpFilename), std::ios::binary | std::ios::trunc);
        if (!out)
            throw FatalError("Couldn't open '" + tmpFilename + "' for output.");
        for (int i : indices)
            out.write(m_points[i].cdata(), m_b.pointSize);
    }
    // Write-then-rename rather than truncating in place: when `key` is a leaf chunk root, its
    // original voxel .bin is still mmap'd by m_points. Renaming over it is safe on POSIX, since the
    // mapping holds the old inode until this ChunkBuilder is destroyed, whereas truncating it
    // would corrupt the points we are reading.
    pdal::FileUtils::renameFile(fullFilename, tmpFilename);
}

// Index a set of prepared chunk plans in parallel, one local octree per chunk. Shared tail of
// both front-ends: runChunkedBuild's leaf-grouped plans and indexChunks's one-file plans.
void indexPlans(const BaseInfo& b, PyramidManager& mgr, std::vector<ChunkPlan>& plans,
    ProgressWriter* progress)
{
    if (plans.empty())
        return;

    // Per-chunk increment: advance the bar from 60% to 90% over the build. Merging finishes it.
    if (progress)
        progress->setIncrement(0.30 / (double)plans.size());

    size_t nthreads = std::thread::hardware_concurrency();
    if (nthreads == 0)
        nthreads = DefaultThreadCount;

    ThreadPool pool(nthreads);
    pool.trap(true);

    for (ChunkPlan& plan : plans)
    {
        ChunkPlan *pp = &plan;
        pool.add([pp, &b, &mgr, progress]()
        {
            // Filenames to delete once the builder has consumed and unmapped them. Exclude the
            // chunk root's own file: if the root is a leaf, ChunkBuilder rewrites that .bin as
            // its output, the subsample the merge reads, so deleting it would drop the output.
            const std::string rootFile = pp->root.toString() + ".bin";
            std::vector<std::string> consumed;
            consumed.reserve(pp->leaves.size());
            for (const auto& l : pp->leaves)
                if (l.second.filename() != rootFile)
                    consumed.push_back(l.second.filename());

            {
                ChunkBuilder cb(b, mgr, std::move(*pp));
                cb.run();
            } // ChunkBuilder and its PointAccessor destroyed here, so leaf files are unmapped.

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
    indexPlans(b, mgr, plans, progress);
}

void indexChunks(const BaseInfo& b, PyramidManager& mgr,
    const std::unordered_map<VoxelKey, FileInfo>& chunkFiles, ProgressWriter* progress)
{
    // Indexing phase: each chunk .bin from Chunking is already exactly one chunk, so build a local
    // octree from each via buildStructure and sampleStructure, like an internal chunk root whose
    // only leaf is itself. root == its own leaf, so writeBin's write-then-rename handles the
    // input==output .bin, and that file is kept for the Merging phase rather than consumed.
    std::vector<ChunkPlan> plans;
    plans.reserve(chunkFiles.size());
    for (const auto& cf : chunkFiles)
    {
        ChunkPlan p;
        p.root = cf.first;
        p.leaves.emplace(cf.first, cf.second);
        plans.push_back(std::move(p));
    }
    std::cerr << "[indexing] chunks=" << plans.size() << "\n";
    indexPlans(b, mgr, plans, progress);
}

} // namespace bu
} // namespace untwine
