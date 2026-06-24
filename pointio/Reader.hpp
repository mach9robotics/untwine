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

#include <pdal/Options.hpp>
#include <pdal/Stage.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/pdal_features.hpp>

#include "untwine/FatalError.hpp"
#include "untwine/FileInfo.hpp"

namespace untwine
{
namespace epf
{

// Build and configure the PDAL reader stage for `fi` (filename + count, plus the LAS-specific
// nosrs/use_eb_vlr/start options). Shared by the EPF and chunker decode passes.
//
// The returned Stage* is owned by `factory`, so the caller must keep `factory` alive for as long as
// the stage is used (e.g. through prepare()/execute()). That's why the factory is passed in by the
// caller rather than created here.
inline pdal::Stage* createReaderStage(pdal::StageFactory& factory, const FileInfo& fi)
{
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

    pdal::Stage *s = factory.createStage(fi.driver);
    if (!s)
        throw FatalError("Couldn't create reader stage for driver '" + fi.driver +
            "' on file '" + fi.filename + "'.");
    s->setOptions(opts);
    return s;
}

} // namespace epf
} // namespace untwine
