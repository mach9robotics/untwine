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

#include <algorithm>

#include <pdal/StageFactory.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/filters/StreamCallbackFilter.hpp>
#include <pdal/pdal_features.hpp>

#include "Distribute.hpp"
#include "ChunkLut.hpp"
#include "Constants.hpp"

#include "../untwine/Common.hpp"
#include "../untwine/FileInfo.hpp"
#include "../untwine/Point.hpp"
#include "../untwine/ProgressWriter.hpp"
#include "../untwine/ThreadPool.hpp"
#include "../untwine/VoxelKey.hpp"
#include "pointio/Cell.hpp"
#include "pointio/EpfTypes.hpp"
#include "pointio/Grid.hpp"
#include "pointio/PointFill.hpp"
#include "pointio/Writer.hpp"

namespace untwine
{
namespace chunker
{

namespace
{

// Decode one file again and route each point through the count grid + LUT into its chunk's .bin.
// Mirrors epf::FileProcessor::run, but the routing key is the chunk root, not the fine voxel.
void distributeFile(FileInfo fi, size_t pointSize, const epf::Grid& grid, const ChunkLut& lut,
    const Transform& xform, epf::Writer *writer, ProgressWriter& progress)
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

    epf::CellMgr cellMgr(pointSize, writer);
    PointCount count = 0;
    epf::Cell *cell = cellMgr.get(VoxelKey());
    epf::PointProcessorPtr ptProcessor = epf::makePointProcessor(fi);

    StreamCallbackFilter f;
    f.setCallback([&cellMgr, &cell, &count, &grid, &lut, &xform, &progress,
                   ptProcessor = ptProcessor.get()](PointRef& point)
    {
        Point p = cell->point();
        ptProcessor->fill(point, p);

        // Quantize identically to the count pass so the cell key matches the LUT.
        p.quantize(xform);

        // Route to the chunk root via the LUT. The count pass saw every point with the same
        // quantization and grid, so the cell is always present. If somehow absent, fall back to
        // the cell itself as its own chunk so no point is ever dropped.
        VoxelKey cellIndex = grid.key(p.x(), p.y(), p.z());
        auto it = lut.cellToRoot.find(cellIndex);
        VoxelKey chunkKey = (it != lut.cellToRoot.end()) ? it->second : cellIndex;

        if (chunkKey != cell->key())
        {
            cell = cellMgr.get(chunkKey, cell);
            cell->copyPoint(p);
        }
        cell->advance();

        count++;
        if (count == ProgressWriter::ChunkSize)
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
        epf::setDimensions(t.layout(), fi);
        f.execute(t);
    }
    catch (const pdal_error& err)
    {
        throw FatalError(err.what());
    }
    progress.update(count);
}

} // unnamed namespace

void distribute(BaseInfo& b, std::vector<FileInfo>& fileInfos, const epf::Grid& grid,
    const ChunkLut& lut, ProgressWriter& progress)
{
    // Largest files first for better thread balance, mirroring EPF.
    std::sort(fileInfos.begin(), fileInfos.end(),
        [](const FileInfo& a, const FileInfo& c) { return a.numPoints > c.numPoints; });

    epf::Writer writer(b.opts.tempDir, epf::NumWriters, b.pointSize);
    // The incrementer is set once by Chunker::run over the whole front-end of 2 * numPoints. Don't
    // reset it here, or the count pass's 0-20% would be lost.

    ThreadPool pool(epf::NumFileProcessors);
    pool.trap(true, "Unknown error in Distribute");
    for (const FileInfo& fi : fileInfos)
    {
        const FileInfo *pfi = &fi;
        pool.add([pfi, &grid, &lut, &progress, &writer, pointSize = b.pointSize, xform = b.xform]()
        {
            distributeFile(*pfi, pointSize, grid, lut, xform, &writer, progress);
        });
    }
    pool.join();
    // stop() blocks until the writer threads finish flushing, then we surface any errors.
    writer.stop();

    StringList errors = pool.clearErrors();
    if (errors.size())
        throw FatalError(errors.front());
}

} // namespace chunker
} // namespace untwine
