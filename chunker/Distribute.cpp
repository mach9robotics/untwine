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
#include <cstring>
#include <vector>

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

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
#include "pointio/Reader.hpp"

namespace untwine
{
namespace chunker
{

namespace
{

// Late materialization: flush a per-file attribute buffer to the shared store with one
// contiguous pwrite. Ids within a file are consecutive, so the destination range
// [firstId * recordSize, ...) is disjoint from every other file's writes.
void flushAttrStore(int attrFd, std::vector<uint8_t>& attrBuf, uint64_t firstId, int recordSize)
{
#ifndef _WIN32
    if (attrBuf.empty())
        return;

    uint64_t off = firstId * (uint64_t)recordSize;
    size_t written = 0;
    while (written < attrBuf.size())
    {
        ssize_t r = ::pwrite(attrFd, attrBuf.data() + written, attrBuf.size() - written,
            (off_t)(off + written));
        if (r <= 0)
            throw FatalError("Failure writing to attribute store.");
        written += r;
    }
    attrBuf.clear();
#endif
}

// Decode one file again and route each point through the count grid + LUT into its chunk's .bin.
// Mirrors epf::FileProcessor::run, but the routing key is the chunk root, not the fine voxel.
//
// attrFd >= 0 selects late materialization: the cell holds only {xyz, rowId} (SlimPointSize
// bytes); the remaining attribute tail is buffered and pwrite'd into the shared attribute store.
void distributeFile(FileInfo fi, size_t pointSize, const epf::Grid& grid, const ChunkLut& lut,
    const Transform& xform, epf::Writer *writer, ProgressWriter& progress, int attrFd,
    int attrRecordSize)
{
    using namespace pdal;

    StageFactory factory;
    Stage *s = epf::createReaderStage(factory, fi);

    epf::CellMgr cellMgr(pointSize, writer);
    PointCount count = 0;
    epf::Cell *cell = cellMgr.get(VoxelKey());
    epf::PointProcessorPtr ptProcessor = epf::makePointProcessor(fi);

    StreamCallbackFilter f;
    if (attrFd >= 0)
    {
        // Late materialization: decode the full point into the wide staging buffer, put
        // {xyz, rowId} in the cell and buffer the attribute tail for the attribute store.
        std::vector<uint8_t> staging(SlimIdOffset + attrRecordSize);
        // Buffer ~2 MB of attribute records; ids within one file are consecutive, so each flush
        // is one contiguous write at a known offset.
        size_t bufRecords = (std::max)((size_t)1, (size_t)(2 * 1024 * 1024) / attrRecordSize);
        size_t attrBufLimit = bufRecords * attrRecordSize;
        std::vector<uint8_t> attrBuf;
        attrBuf.reserve(attrBufLimit);
        uint64_t localRow = 0;
        uint64_t attrFirstId = 0;

        f.setCallback([&cellMgr, &cell, &count, &grid, &lut, &xform, &progress,
                       ptProcessor = ptProcessor.get(), &fi, attrFd, attrRecordSize,
                       &staging, &attrBuf, attrBufLimit, &localRow, &attrFirstId](PointRef& point)
        {
            Point p(staging.data());
            ptProcessor->fill(point, p);

            // Quantize identically to the count pass so the cell key matches the LUT.
            p.quantize(xform);

            // Route to the chunk root via the LUT. The count pass saw every point with the same
            // quantization and grid, so the cell is always present. If somehow absent, fall back
            // to the cell itself as its own chunk so no point is ever dropped.
            VoxelKey cellIndex = grid.key(p.x(), p.y(), p.z());
            auto it = lut.cellToRoot.find(cellIndex);
            VoxelKey chunkKey = (it != lut.cellToRoot.end()) ? it->second : cellIndex;

            if (chunkKey != cell->key())
                cell = cellMgr.get(chunkKey, cell);

            uint64_t id = fi.baseId + localRow++;
            uint8_t *dst = cell->point().data();
            std::memcpy(dst, staging.data(), SlimIdOffset);
            std::memcpy(dst + SlimIdOffset, &id, sizeof(id));
            cell->advance();

            if (attrBuf.empty())
                attrFirstId = id;
            attrBuf.insert(attrBuf.end(), staging.data() + SlimIdOffset,
                staging.data() + SlimIdOffset + attrRecordSize);
            if (attrBuf.size() >= attrBufLimit)
                flushAttrStore(attrFd, attrBuf, attrFirstId, attrRecordSize);

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
        // Write any remaining attribute records.
        flushAttrStore(attrFd, attrBuf, attrFirstId, attrRecordSize);
        progress.update(count);
        return;
    }

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

    int attrFd = -1;
#ifndef _WIN32
    // Late materialization: create the attribute store, sized so that point i's attribute record
    // lives at byte i * attrRecordSize. Each input file owns the disjoint id range
    // [baseId, baseId + numPoints), so the per-file pwrites never overlap.
    if (b.opts.lateMaterialization)
    {
        std::string attrPath = b.opts.tempDir + "/" + AttributeStoreFilename;
        attrFd = ::open(attrPath.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);
        if (attrFd < 0 ||
            ::ftruncate(attrFd, (off_t)(b.numPoints * (uint64_t)b.attrRecordSize)) != 0)
        {
            if (attrFd >= 0)
                ::close(attrFd);
            throw FatalError("Can't create attribute store '" + attrPath + "'.");
        }
    }
#endif

    ThreadPool pool(epf::NumFileProcessors);
    pool.trap(true, "Unknown error in Distribute");
    for (const FileInfo& fi : fileInfos)
    {
        const FileInfo *pfi = &fi;
        pool.add([pfi, &grid, &lut, &progress, &writer, pointSize = b.pointSize, xform = b.xform,
                  attrFd, attrRecordSize = b.attrRecordSize]()
        {
            distributeFile(*pfi, pointSize, grid, lut, xform, &writer, progress, attrFd,
                attrRecordSize);
        });
    }
    pool.join();
#ifndef _WIN32
    // The attribute store is complete once every file processor has flushed.
    if (attrFd >= 0)
        ::close(attrFd);
#endif
    // stop() blocks until the writer threads finish flushing, then we surface any errors.
    writer.stop();

    StringList errors = pool.clearErrors();
    if (errors.size())
        throw FatalError(errors.front());
}

} // namespace chunker
} // namespace untwine
