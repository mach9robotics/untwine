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


#include "untwine/FileInfo.hpp"
#include "Grid.hpp"
#include "Cell.hpp"

namespace untwine
{

class ProgressWriter;

namespace epf
{

class Writer;

// Processes a single input file (FileInfo) and writes data to the Writer.
class FileProcessor
{
public:
    FileProcessor(const FileInfo& fi, size_t pointSize, const Grid& grid, const Transform& xform,
        Writer *writer, ProgressWriter& progress, int attrFd = -1, int attrRecordSize = 0);

    Cell *getCell(const VoxelKey& key);
    void run();

private:
    void flushAttrStore();

    FileInfo m_fi;
    CellMgr m_cellMgr;
    Grid m_grid;
    Transform m_xform;
    ProgressWriter& m_progress;
    // Late materialization (attrFd >= 0): points are decoded into m_staging (full wide
    // layout); {xyz, id} goes to the cell and the attribute tail is buffered in
    // m_attrBuf, flushed with pwrite at (first id * record size).
    int m_attrFd;
    int m_attrRecordSize;
    uint64_t m_localRow {0};
    uint64_t m_attrFirstId {0};
    size_t m_attrBufLimit {0};
    std::vector<uint8_t> m_staging;
    std::vector<uint8_t> m_attrBuf;
};

} // namespace epf
} // namespace untwine
