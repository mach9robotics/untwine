/*****************************************************************************
 *   Copyright (c) 2026, Mach9 Robotics                                      *
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

#include <pdal/util/FileUtils.hpp>

#include "../untwine/Common.hpp"
#include "../untwine/GridKey.hpp"
#include "../untwine/Point.hpp"
#include "../untwine/ProgressWriter.hpp"
#include "../untwine/ThreadPool.hpp"
#include "../untwine/VoxelKey.hpp"

#include "../epf/EpfTypes.hpp"

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

// Leaf cutoff for the in-chunk octree, the same threshold as epf::MaxPointsPerNode, plus a depth
// guard so pathological coincident points can't recurse forever.
constexpr size_t LeafPointBudget = epf::MaxPointsPerNode;
constexpr int MaxBuildLevel = 30;

// Fallback worker-thread count when std::thread::hardware_concurrency() reports 0.
constexpr unsigned DefaultThreadCount = 8;

// One chunk to build: an internal subtree root plus its descendant leaf files.
struct ChunkPlan
{
    VoxelKey root;
    std::unordered_map<VoxelKey, FileInfo> leaves;
};

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

    // Build a self-contained chunk for `indices` and hand it to the ChunkWriter for deferred
    // compression and writing.
    void emit(const VoxelKey& key, const std::vector<int>& indices);
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

    // Build a self-contained chunk: a shared PointTable plus a view of copied point data. The view
    // owns its data, so compression and writing run on a ChunkWriter thread after this ChunkBuilder
    // and its mmap'd point files are gone. Stats accumulation and logOctant happen there too.
    auto table = std::make_shared<PointTable>();
    IndexedStats stats;
    DimTypeList extraDims;
    DimInfoList dims = m_b.dimInfo;
    for (FileDimInfo& fdi : dims)
    {
        fdi.dim = table->layout()->registerOrAssignDim(fdi.name, fdi.type);
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
    table->finalize();

    PointViewPtr view(new PointView(*table));
    PointId pointId = 0;
    for (int idx : indices)
    {
        char *base = m_points[idx].cdata();
        for (const FileDimInfo& fdi : dims)
            view->setField(fdi.dim, fdi.type, pointId,
                reinterpret_cast<void *>(base + fdi.offset));
        pointId++;
    }

    m_mgr.enqueueChunk({key, table, view, extraDims, std::move(stats), indices.size()});
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

// Index a set of prepared chunk plans in parallel, one local octree per chunk. The tail of the
// chunker indexing phase: indexChunks builds one plan per chunk .bin and calls this.
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
