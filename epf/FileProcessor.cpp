/*****************************************************************************
 *   Copyright (c) 2020, Hobu, Inc. (info@hobu.co)                           *
 *                                                                           *
 *   All rights reserved.                                                    *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 ****************************************************************************/


#include "FileProcessor.hpp"
#include "../untwine/ProgressWriter.hpp"

#ifndef _WIN32
#include <unistd.h>
#endif

#include <pdal/pdal_features.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/filters/StreamCallbackFilter.hpp>

namespace untwine
{
namespace epf
{

namespace
{

// The dimension IDs need to come from the *source* layout because it's possible in the
// case of user-defined dimensions that the IDs could vary for the same-named dimension
// in different input files. The "dim" field represents the ID of the dimension
// we're reading from. "offset" is the corresponding notion in the output packed point data.
void setDimensions(pdal::PointLayoutPtr layout, FileInfo& fi)
{
    for (FileDimInfo& di : fi.dimInfo)
    {
        di.dim = layout->findDim(di.name);
        assert(di.dim != pdal::Dimension::Id::Unknown);

        // If we have a bit offset, then we're really writing to the classflags field in the
        // output point. Fetch that offset and set it. We will probably do this several
        // times (once for each bit), but the value should always be the same.
        if (di.shift != -1)
            fi.untwineBitsOffset = di.offset;
    }
}

} // unnamed namespace


FileProcessor::FileProcessor(const FileInfo& fi, size_t pointSize, const Grid& grid,
        const Transform& xform, Writer *writer, ProgressWriter& progress, int attrFd,
        int attrRecordSize) :
    m_fi(fi), m_cellMgr(pointSize, writer), m_grid(grid), m_xform(xform), m_progress(progress),
    m_attrFd(attrFd), m_attrRecordSize(attrRecordSize)
{
    if (m_attrFd >= 0)
    {
        m_staging.resize(SlimIdOffset + m_attrRecordSize);
        // Buffer ~2 MB of attribute records; ids within a FileProcessor are consecutive,
        // so each flush is one contiguous write at a known offset.
        size_t bufRecords = (std::max)((size_t)1, (size_t)(2 * 1024 * 1024) / m_attrRecordSize);
        m_attrBufLimit = bufRecords * m_attrRecordSize;
        m_attrBuf.reserve(m_attrBufLimit);
    }
}

void FileProcessor::flushAttrStore()
{
#ifndef _WIN32
    if (m_attrBuf.empty())
        return;

    uint64_t off = m_attrFirstId * (uint64_t)m_attrRecordSize;
    size_t written = 0;
    while (written < m_attrBuf.size())
    {
        ssize_t r = ::pwrite(m_attrFd, m_attrBuf.data() + written, m_attrBuf.size() - written,
            (off_t)(off + written));
        if (r <= 0)
            throw FatalError("Failure writing to attribute store.");
        written += r;
    }
    m_attrBuf.clear();
#endif
}

class BasePointProcessor
{
public:
    BasePointProcessor(const FileInfo& fi) : m_fi(fi)
    {}
    virtual ~BasePointProcessor()
    {}

    virtual void fill(const pdal::PointRef& src, Point& dst) = 0;

protected:
    const FileInfo& m_fi;
};
using PointProcessorPtr = std::unique_ptr<BasePointProcessor>;

// These processors could probably be improved performance-wise by breaking the dimensions
// up into types in the ctor to avoid the conditionals in fill().
// Could also make FileProcessor take these as a template type to avoid having fill() be
// virtual.
class StdPointProcessor : public BasePointProcessor
{
public:
    using BasePointProcessor::BasePointProcessor;

    void fill(const pdal::PointRef& src, Point& dst) override
    {
        uint8_t untwineBits = 0;
        for (const FileDimInfo& fdi : m_fi.dimInfo)
        {
            if (fdi.shift == -1)
                src.getField(reinterpret_cast<char *>(dst.data() + fdi.offset),
                    fdi.dim, fdi.type);
            else
                untwineBits |= (src.getFieldAs<uint8_t>(fdi.dim) << fdi.shift);
        }

        // We pack all the bitfields into the "untwine bits" field.
        if (m_fi.untwineBitsOffset > -1)
            memcpy(dst.data() + m_fi.untwineBitsOffset, &untwineBits, 1);
    }
};

class LegacyLasPointProcessor : public BasePointProcessor
{
public:
    using BasePointProcessor::BasePointProcessor;

    void fill(const pdal::PointRef& src, Point& dst) override
    {
        uint8_t untwineBits = 0;
        for (const FileDimInfo& fdi : m_fi.dimInfo)
        {
            if (fdi.dim == pdal::Dimension::Id::Classification)
            {
                uint8_t classification = src.getFieldAs<uint8_t>(fdi.dim);
                if (classification == 12)
                    untwineBits |= 0x08;  // Set the overlap bit.
                untwineBits |= (classification >> 5);
                classification &= 0x1F;
                memcpy(dst.data() + fdi.offset, &classification, 1);
            }
            else if (fdi.shift == -1)
                src.getField(reinterpret_cast<char *>(dst.data() + fdi.offset),
                    fdi.dim, fdi.type);
            else
                untwineBits |= (src.getFieldAs<uint8_t>(fdi.dim) << fdi.shift);
        }

        // We pack all the bitfields into the "untwine bits" field.
        if (m_fi.untwineBitsOffset > -1)
            memcpy(dst.data() + m_fi.untwineBitsOffset, &untwineBits, 1);
    }
};


void FileProcessor::run()
{
    pdal::Options opts;
    opts.add("filename", m_fi.filename);
    opts.add("count", m_fi.numPoints);
    if (m_fi.driver == "readers.las")
    {
        opts.add("nosrs", m_fi.no_srs);
        opts.add("use_eb_vlr", "true");
#ifdef PDAL_LAS_START
        opts.add("start", m_fi.start);
#endif
    }

    pdal::StageFactory factory;
    pdal::Stage *s = factory.createStage(m_fi.driver);
    s->setOptions(opts);

    PointCount count = 0;

    // We need to move the data from the PointRef to some output buffer. We copy the data
    // to the end of the *last* output buffer we used in hopes that it's the right one.
    // If it's not we lose and we're forced to move that data to the another cell,
    // which then becomes the active cell.

    // This is some random cell that ultimately won't get used, but it contains a buffer
    // into which we can write data.
    Cell *cell = m_cellMgr.get(VoxelKey());

    PointProcessorPtr ptProcessor;
    if (m_fi.driver == "readers.las" && m_fi.fileVersion < 14)
        ptProcessor = std::make_unique<LegacyLasPointProcessor>(m_fi);
    else
        ptProcessor = std::make_unique<StdPointProcessor>(m_fi);

    pdal::StreamCallbackFilter f;
    if (m_attrFd >= 0)
    {
        // Late materialization: decode the full point into the wide staging buffer, put
        // {xyz, rowId} in the cell and buffer the attribute tail for the attribute store.
        // Nothing in-progress ever lives in a cell buffer, so no flush-exclusion is needed.
        f.setCallback([this, &count, &cell, ptProcessor = ptProcessor.get()](pdal::PointRef& point)
            {
                Point p(m_staging.data());
                ptProcessor->fill(point, p);

                // Quantize the point value to what will be stored in the output file
                // so that we're sure that the WRITTEN point value is in the expected
                // location in the grid.
                p.quantize(m_xform);

                VoxelKey cellIndex = m_grid.key(p.x(), p.y(), p.z());
                if (cellIndex != cell->key())
                    cell = m_cellMgr.get(cellIndex);

                uint64_t id = m_fi.baseId + m_localRow++;
                uint8_t *dst = cell->point().data();
                memcpy(dst, m_staging.data(), SlimIdOffset);
                memcpy(dst + SlimIdOffset, &id, sizeof(id));
                cell->advance();

                if (m_attrBuf.empty())
                    m_attrFirstId = id;
                m_attrBuf.insert(m_attrBuf.end(), m_staging.data() + SlimIdOffset,
                    m_staging.data() + SlimIdOffset + m_attrRecordSize);
                if (m_attrBuf.size() >= m_attrBufLimit)
                    flushAttrStore();

                count++;
                if (count == ProgressWriter::ChunkSize)
                {
                    m_progress.update();
                    count = 0;
                }

                return true;
            }
        );
    }
    else
    f.setCallback([this, &count, &cell, ptProcessor = ptProcessor.get()](pdal::PointRef& point)
        {
            // Write the data into the point buffer in the cell.  This is the *last*
            // cell buffer that we used. We're hoping that it's the right one.
            Point p = cell->point();
            ptProcessor->fill(point, p);

            // Quantize the point value to what will be stored in the output file
            // so that we're sure that the WRITTEN point value is in the expected
            // location in the grid.
            p.quantize(m_xform);

            // Find the actual cell that this point belongs in. If it's not the one
            // we chose, copy the data to the correct cell.
            VoxelKey cellIndex = m_grid.key(p.x(), p.y(), p.z());
            if (cellIndex != cell->key())
            {
                // Make sure that we exclude the current cell from any potential flush so
                // that no other thread can claim its data buffer and overwrite it before
                // we have a chance to copy from it in copyPoint().
                cell = m_cellMgr.get(cellIndex, cell);
                cell->copyPoint(p);
            }
            // Advance the cell - move the buffer pointer so when we refer to the cell's
            // point, we're referring to the next location in the cell's buffer.
            cell->advance();
            count++;

            if (count == ProgressWriter::ChunkSize)
            {
                m_progress.update();
                count = 0;
            }

            return true;
        }
    );
    f.setInput(*s);

    pdal::FixedPointTable t(1000);

    try
    {
        f.prepare(t);
        setDimensions(t.layout(), m_fi);
        f.execute(t);
    }
    catch (const pdal::pdal_error& err)
    {
        throw FatalError(err.what());
    }

    // Late materialization: write any remaining attribute records.
    flushAttrStore();

    // We normally call update for every CountIncrement points, but at the end, just
    // tell the progress writer the number that we've done since the last update.
    m_progress.update(count);
}

} // namespace epf
} // namespace untwine
