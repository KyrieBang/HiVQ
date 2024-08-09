#ifndef CORE_CREATOR_HPP_
#define CORE_CREATOR_HPP_

#include "core.hpp"
#include "core_iterator.hpp"
#include "geo.hpp"
#include "mem.hpp"
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

namespace HiGIS::Core {

class PolygonCreator;
class LineStringCreator;
class MultiPointCreator;
class PointCreator;

void PushExtent(MemPool<Extent> &out, double *points, uint64_t n);
void FlushContent(GeoBuffer &buffer, PolygonCreator &c);
void FlushContent(GeoBuffer &buffer, LineStringCreator &c);
void FlushContent(GeoBuffer &buffer, MultiPointCreator &c);
/**
 * \brief Flush content of Point to data.
 * \ingroup core_creator
 */
void Flush(GeoBuffer &buffer, GeoSummary &meta, PointCreator &c);

/**
 * \brief Flush content of Polygon, LineString or MultiPoint to data.
 * \ingroup core_creator
 */
template <typename CreatorT>
void Flush(GeoBuffer &buffer, GeoSummary &summary, Mix<FeatureExtents> &stat, CreatorT &c) {
    buffer.Push(c.id == INFINITY ? summary.feature_count : c.id); // id
    if (c.NodeData().Size() == 0) {
        c.FlushEmpty(buffer);
    } else {
        FlushContent(buffer, c);
    }
    buffer.Push(TERMINATOR);
    summary.feature_count += 1;
    PushExtent(stat.extents, c.NodeData().Data(), c.NodeData().Size());
    if (c.NodeData().Size() > 0) {
        UpdateExtent(summary.extent, stat.extents[-1]);
    }
}

template <typename TIter> void CopyFeature(GeoBuffer &out, TIter iter) {
    out.Copy(iter.Ptr().begin, iter.Length());
}

/** \brief Create polygon and flush to pool */
class PolygonCreator {
 public:
    uint64_t id = INFINITY; // infinity means auto

    struct Options {
        bool auto_close_ring = true;
        bool auto_clock_direction = true;
    };
    explicit PolygonCreator(Options options);
    PolygonCreator() : PolygonCreator(Options()) {}
    MemPool<uint64_t> Index2() { return index2_; }    // for parts
    MemPool<uint64_t> Index1() { return index1_; }    // for rings
    MemPool<double> NodeData() { return node_data_; } // for nodes
    PolygonCreator &AddPart();
    PolygonCreator &AddRing();
    PolygonCreator &AddPoint(double x, double y);
    PolygonCreator &FixRing();
    PolygonAccessor Parse();
    void Clear();
    // syntax sugar
    PolygonCreator &Part() { return AddPart(); }
    PolygonCreator &Ring() { return AddRing(); }
    PolygonCreator &operator()(double x, double y) { return AddPoint(x, y); }
    void Flush(GeoData &out) { Core::Flush(out, out, out, *this); }
    static void FlushEmpty(GeoBuffer &out);

 private:
    MemPool<uint64_t> index2_;  // for parts
    MemPool<uint64_t> index1_;  // for rings
    MemPool<double> node_data_; // for nodes
    Options options_;
};

/** \brief Create linestring and flush to pool */
class LineStringCreator {
 public:
    uint64_t id = INFINITY; // infinity means auto
    LineStringCreator() { index1_.Push(0); }
    MemPool<uint64_t> Index1() { return index1_; }    // for lines
    MemPool<double> NodeData() { return node_data_; } // for nodes
    LineStringCreator &AddLine();
    LineStringCreator &AddPoint(double x, double y);
    LineStringAccessor Parse();
    void Clear();
    // syntax sugar
    LineStringCreator &Line() { return AddLine(); }
    LineStringCreator &operator()(double x, double y) { return AddPoint(x, y); }
    void Flush(GeoData &out) { Core::Flush(out, out, out, *this); }
    static void FlushEmpty(GeoBuffer &out);

 private:
    MemPool<uint64_t> index1_;  // for rings
    MemPool<double> node_data_; // for nodes
};

/** \brief Create multipoint and flush to pool */
class MultiPointCreator {
 public:
    uint64_t id = INFINITY;                           // infinity means auto
    MemPool<double> NodeData() { return node_data_; } // for nodes
    MultiPointCreator &AddPoint(double x, double y);
    MultiPointAccessor Parse();
    void Clear();
    // syntax sugar
    MultiPointCreator &operator()(double x, double y) { return AddPoint(x, y); }
    void Flush(GeoData &out) { Core::Flush(out, out, out, *this); }
    static void FlushEmpty(GeoBuffer &out);

 private:
    MemPool<double> node_data_; // for nodes
};

/** \brief Create point and flush to pool */
class PointCreator {
 public:
    uint64_t id = INFINITY; // infinity means auto
    PointCreator(double x, double y) : x(x), y(y) {}
    /** \brief Constructor for empty point, always pass false to get an empty point */
    PointCreator(bool has_value) : empty_(!has_value) {}
    PointAccessor Parse();
    void Flush(GeoData &out) { Core::Flush(out, out, *this); }
    void Clear();
    bool Empty() { return empty_; }
    static void FlushEmpty(GeoBuffer &out);
    double x;
    double y;

 private:
    bool empty_ = false;
};

void FlushEmpty(GeoBuffer &out, GeoType type);

} // namespace HiGIS::Core

#endif // CORE_CREATOR_HPP_

/**
 * \defgroup core_creator Core Creator
 * \brief Provides classes to create new features, much easier than directly write binary memory
 *
 * Creators are append-only facilities to gradually append components of a feature, and then flush
 * into a GeoData or components of geo datasets. The API style of creators are fluent chains,
 * providing syntax sugars to be more succinct.
 *
 * Creators can generate Iterators (core_iterator) to access the half-baked geometry.
 *
 * Example (PolygonCreator)
 *
 * \copydoc PolygonCreator
 */

/**
 * \class HiGIS::Core::PolygonCreator
 * \ingroup core_creator
 *
 * \snippet test_creator.cpp PolygonCreator
 *
 * More succinct APIs for simpler geometries
 * \snippet test_creator.cpp PolygonCreatorSimple
 *
 * Flush in a more flexible manner
 * \snippet test_creator.cpp PolygonCreatorFlush2
 *
 * Access the half-baked geometry:
 * \snippet test_creator.cpp PolygonCreatorAccess
 */

/**
 * \class HiGIS::Core::LineStringCreator
 * \ingroup core_creator
 *
 * \snippet test_creator.cpp LineStringCreator
 *
 * Access the half baked geometry
 * \snippet test_creator.cpp LineStringCreatorAccess
 */

/**
 * \class HiGIS::Core::MultiPointCreator
 * \ingroup core_creator
 *
 * \snippet test_creator.cpp MultiPointCreator
 *
 * Access the half baked geometry
 * \snippet test_creator.cpp MultiPointCreatorAccess
 */

/**
 * \class HiGIS::Core::PointCreator
 * \ingroup core_creator
 *
 * \snippet test_creator.cpp PointCreator
 *
 * Access the half baked geometry
 * \snippet test_creator.cpp PointCreatorAccess
 *
 */