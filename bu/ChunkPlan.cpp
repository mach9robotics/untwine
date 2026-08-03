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

#include <functional>

#include "ChunkPlan.hpp"

namespace untwine
{
namespace bu
{

std::vector<VoxelKey> planChunkRoots(const std::unordered_map<VoxelKey, uint64_t>& leafCounts,
    uint64_t target)
{
    const VoxelKey root;

    // Total leaf points in the subtree under every key (the leaf itself and each ancestor).
    std::unordered_map<VoxelKey, uint64_t> cnt;
    for (const auto& p : leafCounts)
    {
        VoxelKey k = p.first;
        while (true)
        {
            cnt[k] += p.second;
            if (k == root)
                break;
            k = k.parent();
        }
    }

    std::vector<VoxelKey> roots;
    std::function<void(const VoxelKey&)> descend = [&](const VoxelKey& k)
    {
        auto ci = cnt.find(k);
        if (ci == cnt.end() || ci->second == 0)
            return;
        bool isLeaf = leafCounts.count(k) > 0;
        bool fits = ci->second <= target;
        if (k != root && (isLeaf || fits))
        {
            // Internal nodes that fit the budget become chunk roots. A leaf that exceeds the
            // budget also becomes a buildable root, since ChunkBuilder splits it in RAM and there
            // is no longer an on-disk reprocessor to pre-split it. A leaf that fits is left
            // standalone for the merge.
            if (!isLeaf || ci->second > target)
                roots.push_back(k);
            return;
        }
        for (int i = 0; i < 8; ++i)
            descend(k.child(i));
    };
    descend(root);
    return roots;
}

int localGridDepth(uint64_t numPoints, uint64_t leafBudget, int maxDepthBelowRoot)
{
    if (leafBudget == 0)
        leafBudget = 1;

    // Smallest d with 8^d >= ceil(numPoints / leafBudget): a uniform spread then puts at most
    // leafBudget points per finest cell. One extra level of margin so mildly clustered cells fit.
    uint64_t ratio = (numPoints + leafBudget - 1) / leafBudget;
    int d = 1;
    while (d < 20 && (uint64_t(1) << (3 * d)) < ratio)
        ++d;
    ++d;  // margin

    if (maxDepthBelowRoot >= 1 && d > maxDepthBelowRoot)
        d = maxDepthBelowRoot;
    if (d < 1)
        d = 1;
    return d;
}

} // namespace bu
} // namespace untwine
