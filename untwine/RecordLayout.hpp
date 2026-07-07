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

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pdal/Dimension.hpp>

#include "FileDimInfo.hpp"

namespace untwine
{

// A dimension's location within a packed temp record: its byte offset and stored type.
//
// A Field is only produced by RecordLayout, which validates offset + size(type) <= pointSize when
// it is built, so a Field handed to PackedPoint::get<T> is guaranteed in-bounds. offset < 0 marks
// an absent dimension; reads through it yield T{}, mirroring PDAL's getFieldAs on a missing dim.
struct Field
{
    int offset {-1};
    pdal::Dimension::Type type {pdal::Dimension::Type::None};

    bool present() const { return offset >= 0; }
};

namespace detail
{

// Read sizeof(Src) bytes at p as Src, then numerically cast to T. memcpy (not a pointer cast)
// handles unaligned offsets and avoids a strict-aliasing violation. Compiles to a load + cast.
template<typename Src, typename T>
inline T readAs(const char* p)
{
    Src v;
    std::memcpy(&v, p, sizeof(v));
    return static_cast<T>(v);
}

} // namespace detail

class RecordLayout;

// A read-only view of one packed point record, bound to its RecordLayout. It absorbs the offset
// arithmetic and runtime type dispatch that would otherwise be hand-threaded at every call site:
// get<T>(field) reads the value stored at field.offset (interpreted as field.type) and numerically
// casts it to T. All field access to the packed .bin record format should go through this rather
// than reading rec + offset directly, so a wrong offset/type can't leak in silently.
class PackedPoint
{
public:
    PackedPoint(const char* rec, const RecordLayout& layout) : m_rec(rec), m_layout(&layout)
    {}

    // Read `f` from this record and numerically cast to T. Absent field -> T{}. `f` must come from
    // this record's layout; its bounds were validated when the layout was built, so the read below
    // is in-bounds by construction (no per-access bounds check needed).
    template<typename T>
    T get(const Field& f) const
    {
        using Type = pdal::Dimension::Type;

        if (!f.present())
            return T{};

        // The stored type is only known at runtime, and converting its bytes to T is a numeric
        // cast (not a reinterpret), so we dispatch on it. This is getFieldAs<T> minus the layout
        // lookup — reachable only with a validated Field, never a hand-passed offset/type pair.
        const char* p = m_rec + f.offset;
        switch (f.type)
        {
        case Type::Unsigned8:  return detail::readAs<uint8_t,  T>(p);
        case Type::Unsigned16: return detail::readAs<uint16_t, T>(p);
        case Type::Unsigned32: return detail::readAs<uint32_t, T>(p);
        case Type::Unsigned64: return detail::readAs<uint64_t, T>(p);
        case Type::Signed8:    return detail::readAs<int8_t,   T>(p);
        case Type::Signed16:   return detail::readAs<int16_t,  T>(p);
        case Type::Signed32:   return detail::readAs<int32_t,  T>(p);
        case Type::Signed64:   return detail::readAs<int64_t,  T>(p);
        case Type::Float:      return detail::readAs<float,    T>(p);
        case Type::Double:     return detail::readAs<double,   T>(p);
        case Type::None:       return T{};
        }
        return T{};
    }

    // Convenience: resolve `id` -> Field via the layout, then read. Does a hash lookup per call, so
    // hot per-point loops should resolve a Field once (RecordLayout::field) and pass it to
    // get<T>(Field); this overload is for cold or one-off reads. Defined below RecordLayout.
    template<typename T>
    T get(pdal::Dimension::Id id) const;

    // Copy `f`'s raw stored bytes (size(f.type) bytes) verbatim into dst — no numeric cast. For the
    // LAZ extra-byte passthrough, which stores each value bit-for-bit. No-op if `f` is absent, so
    // the caller must zero dst first to get the absent-dim default.
    void copyRaw(const Field& f, void* dst) const
    {
        if (f.present())
            std::memcpy(dst, m_rec + f.offset, pdal::Dimension::size(f.type));
    }

    const char* data() const { return m_rec; }

private:
    const char* m_rec;
    const RecordLayout* m_layout;
};

// Owns the packed-record byte layout for one job: the record size (pointSize) plus the
// id/name -> Field table, built once from a DimInfoList. The construction validates the layout
// (every field within pointSize, no two fields overlapping), turning "trust the offset math" into
// "proved it once." It then hands out Fields (resolved locations) and PackedPoint views;
// point(base, i) absorbs the i * pointSize stride into the packed buffer.
//
// Any reader of the .bin format should build one of these from the job's dimInfo (BaseInfo) and
// route all field access through the PackedPoint views it produces.
class RecordLayout
{
public:
    RecordLayout(const DimInfoList& dims, size_t pointSize) : m_pointSize(pointSize)
    {
        for (const FileDimInfo& fdi : dims)
        {
            Field f{fdi.offset, fdi.type};
            assert((!f.present() ||
                    f.offset + (int)pdal::Dimension::size(f.type) <= (int)m_pointSize) &&
                "packed record field runs past the end of the record");
            // Std fields are looked up by canonical id; the packed "untwine bits" field has a
            // dynamically-assigned id, so it is also keyed by name (see ChunkWriter::resolveStdLocs).
            m_byId[fdi.dim] = f;
            m_byName[fdi.name] = f;
        }
        assertNoOverlaps(dims);
    }

    size_t pointSize() const { return m_pointSize; }

    // Resolved location for `id` (absent Field if the dimension isn't stored). Resolve once and
    // reuse across a point loop to keep the read out of the hash table.
    Field field(pdal::Dimension::Id id) const
    {
        auto it = m_byId.find(id);
        return it != m_byId.end() ? it->second : Field{};
    }

    // Resolved location by dimension name (for the packed "untwine bits" field, whose id is not
    // canonical).
    Field field(const std::string& name) const
    {
        auto it = m_byName.find(name);
        return it != m_byName.end() ? it->second : Field{};
    }

    // A view of point `i` within the packed buffer `base` (base == records.data()). This is the
    // single place the i * pointSize stride lives.
    PackedPoint point(const char* base, size_t i) const
    {
        return PackedPoint(base + i * m_pointSize, *this);
    }

private:
    // Debug-only: assert no two stored fields' byte ranges overlap. The packed layout comes from a
    // PDAL PointLayout (prep/FilePrep::fillMetadata), where offsets are contiguous by construction,
    // so this only guards against a future producer that builds the layout by hand.
    void assertNoOverlaps(const DimInfoList& dims) const
    {
#ifndef NDEBUG
        std::vector<std::pair<int, int>> ranges;   // [offset, offset + size)
        for (const FileDimInfo& fdi : dims)
            if (fdi.offset >= 0 && fdi.type != pdal::Dimension::Type::None)
                ranges.emplace_back(fdi.offset, fdi.offset + (int)pdal::Dimension::size(fdi.type));
        std::sort(ranges.begin(), ranges.end());
        for (size_t i = 1; i < ranges.size(); ++i)
            assert(ranges[i].first >= ranges[i - 1].second &&
                "packed record fields overlap");
#else
        (void)dims;
#endif
    }

    size_t m_pointSize;
    std::unordered_map<pdal::Dimension::Id, Field> m_byId;
    std::unordered_map<std::string, Field> m_byName;
};

template<typename T>
inline T PackedPoint::get(pdal::Dimension::Id id) const
{
    return get<T>(m_layout->field(id));
}

} // namespace untwine
