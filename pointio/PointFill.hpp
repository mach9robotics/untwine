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

#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <memory>

#include <pdal/PointLayout.hpp>
#include <pdal/PointRef.hpp>

#include "../untwine/FileInfo.hpp"
#include "../untwine/Point.hpp"

// Serializer: PDAL point -> untwine's packed fixed-size .bin record. Shared between the EPF
// FileProcessor, which bins into per-voxel .bins, and the chunker Distribute pass, which routes
// into per-chunk .bins. The byte layout is defined by the BaseInfo dim layout: the FileDimInfo
// offsets plus the packed "untwine bits" field. Any producer of the .bin format must use this fill.

namespace untwine
{
namespace epf
{

// Resolve each of `fi`'s dimensions against the source layout. The dimension IDs need to come from
// the source layout because for user-defined dimensions the IDs could vary for the same-named
// dimension across input files. The "dim" field is the ID we read from; "offset" is the
// corresponding location in the output packed point.
//
// layout: the source point layout to resolve dimension IDs against.
// fi:     the file info whose dimInfo and untwineBitsOffset are filled in place.
inline void setDimensions(pdal::PointLayoutPtr layout, FileInfo& fi)
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

// Abstract base for the per-file point serializers. A processor copies the dimensions named in its
// FileInfo out of a source PDAL point into a packed output record, using the file's resolved
// dim->offset layout. Subclasses differ only in how they pack the LAS classification bit fields.
class BasePointProcessor
{
public:
    BasePointProcessor(const FileInfo& fi) : m_fi(fi)
    {}
    virtual ~BasePointProcessor()
    {}

    // Fill packed output record `dst` from source point `src`, per m_fi's dimension layout.
    virtual void fill(const pdal::PointRef& src, Point& dst) = 0;

protected:
    const FileInfo& m_fi;
};
using PointProcessorPtr = std::unique_ptr<BasePointProcessor>;

// Serializer for LAS 1.4+ and non-LAS files: copies each dimension through to its offset, OR-ing
// any bit-field dimensions into the packed "untwine bits" byte.
// These processors could probably be improved performance-wise by breaking the dimensions
// up into types in the ctor to avoid the conditionals in fill().
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

// Serializer for pre-1.4 LAS, where classification also carries the overlap flag (class 12) and
// the high classification bits, which are split out into the packed "untwine bits" byte here.
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

// Pick the fill processor for a file. Pre-1.4 LAS packs the classification overlap and synthetic
// bits differently, so the choice depends on the file version.
//
// fi:      the file info, read for its driver and file version.
// Returns: a processor that fills packed records for this file's layout.
inline PointProcessorPtr makePointProcessor(const FileInfo& fi)
{
    if (fi.driver == "readers.las" && fi.fileVersion < 14)
        return std::make_unique<LegacyLasPointProcessor>(fi);
    return std::make_unique<StdPointProcessor>(fi);
}

} // namespace epf
} // namespace untwine
