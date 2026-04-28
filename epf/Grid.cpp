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


#include <cmath>
#include <cstdint>

#include "Epf.hpp"
#include "EpfTypes.hpp"
#include "Grid.hpp"

using namespace pdal;

namespace untwine
{
namespace epf
{

int Grid::calcLevel(uint64_t numPoints, bool cubic, double fillRatio) const
{
    double xside = m_bounds.maxx - m_bounds.minx;
    double yside = m_bounds.maxy - m_bounds.miny;
    double zside = m_bounds.maxz - m_bounds.minz;

    double side = (std::max)(xside, (std::max)(yside, zside));

    // Adjust point count to reflect actual density within the occupied volume.
    // If data only fills a fraction of the cube, points are denser in the occupied
    // region than a uniform distribution would suggest. This avoids underestimating
    // the grid level for disjoint or spatially sparse inputs.
    double pts = (double)numPoints;
    if (fillRatio > 0 && fillRatio < 1.0)
        pts = numPoints / fillRatio;

    int level = 0;
    while (pts > MaxPointsPerNode)
    {
        if (cubic)
        {
            if (xside >= side)
                pts /= 2;
            if (yside >= side)
                pts /= 2;
            if (zside >= side)
                pts /= 2;
        }
        else
            pts /= 8;
        side /= 2;
        level++;
    }

    return level;
}

void Grid::resetLevel(int level)
{
    // We have to have at least level 1 or things break when sampling.
    m_maxLevel = (std::max)(level, 1);
    m_gridSize = (int)std::pow(2, m_maxLevel);
    m_xsize = (m_bounds.maxx - m_bounds.minx) / m_gridSize;
    m_ysize = (m_bounds.maxy - m_bounds.miny) / m_gridSize;
    m_zsize = (m_bounds.maxz - m_bounds.minz) / m_gridSize;
}

VoxelKey Grid::key(double x, double y, double z) const
{
    int xi = (int)std::floor((x - m_bounds.minx) / m_xsize);
    int yi = (int)std::floor((y - m_bounds.miny) / m_ysize);
    int zi = (int)std::floor((z - m_bounds.minz) / m_zsize);
    xi = (std::min)((std::max)(0, xi), m_gridSize - 1);
    yi = (std::min)((std::max)(0, yi), m_gridSize - 1);
    zi = (std::min)((std::max)(0, zi), m_gridSize - 1);

    return VoxelKey(xi, yi, zi, m_maxLevel);
}

} // namespace epf
} // namespace untwine
