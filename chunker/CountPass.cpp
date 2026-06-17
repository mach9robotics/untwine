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

#include <mutex>

#include <pdal/StageFactory.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/filters/StreamCallbackFilter.hpp>
#include <pdal/pdal_features.hpp>

#include "CountPass.hpp"
#include "Constants.hpp"

#include "../untwine/Common.hpp"
#include "../untwine/FileInfo.hpp"
#include "../untwine/Point.hpp"
#include "../untwine/ProgressWriter.hpp"
#include "../untwine/ThreadPool.hpp"
#include "../epf/Grid.hpp"

namespace untwine
{
namespace chunker
{

namespace
{

// Decode one file and tally per-cell counts. Only X/Y/Z are read. Quantization mirrors the
// distribute pass via Point::quantize against the output transform so the cell keys match.
std::unordered_map<VoxelKey, uint64_t> countFile(const FileInfo& fi, const epf::Grid& grid,
    const Transform& xform, ProgressWriter& progress)
{
    using namespace pdal;

    pdal::Options opts;
    opts.add("filename", fi.filename);
    opts.add("count", fi.numPoints);
    if (fi.driver == "readers.las")
    {
        opts.add("nosrs", fi.no_srs);
        opts.add("use_eb_vlr", "true");
#ifdef PDAL_LAS_START
        opts.add("start", fi.start);
#endif
    }

    StageFactory factory;
    Stage *s = factory.createStage(fi.driver);
    s->setOptions(opts);

    std::unordered_map<VoxelKey, uint64_t> counts;
    PointCount count = 0;
    StreamCallbackFilter f;
    f.setCallback([&counts, &grid, &xform, &count, &progress](PointRef& point)
    {
        double xyz[3] = {
            point.getFieldAs<double>(Dimension::Id::X),
            point.getFieldAs<double>(Dimension::Id::Y),
            point.getFieldAs<double>(Dimension::Id::Z)
        };
        Point p(reinterpret_cast<char *>(xyz));
        p.quantize(xform);
        counts[grid.key(p.x(), p.y(), p.z())]++;

        if (++count == ProgressWriter::ChunkSize)
        {
            progress.update();
            count = 0;
        }
        return true;
    });
    f.setInput(*s);

    FixedPointTable t(StreamPointBatch);
    try
    {
        f.prepare(t);
        f.execute(t);
    }
    catch (const pdal_error& err)
    {
        throw FatalError(err.what());
    }
    progress.update(count);
    return counts;
}

} // unnamed namespace

std::unordered_map<VoxelKey, uint64_t> countPoints(std::vector<FileInfo>& fileInfos,
    const epf::Grid& grid, const Transform& xform, std::size_t numThreads,
    ProgressWriter& progress)
{
    std::unordered_map<VoxelKey, uint64_t> total;
    std::mutex mtx;

    ThreadPool pool(numThreads);
    pool.trap(true);
    for (FileInfo& fi : fileInfos)
    {
        FileInfo *pfi = &fi;
        pool.add([pfi, &grid, &xform, &total, &mtx, &progress]()
        {
            std::unordered_map<VoxelKey, uint64_t> local = countFile(*pfi, grid, xform, progress);
            std::lock_guard<std::mutex> lock(mtx);
            for (const auto& c : local)
                total[c.first] += c.second;
        });
    }
    pool.await();
    if (pool.hasErrors())
    {
        std::vector<std::string> errs = pool.clearErrors();
        throw FatalError(errs.empty() ? std::string("Count pass failed.") : errs.front());
    }
    pool.join();
    return total;
}

} // namespace chunker
} // namespace untwine
