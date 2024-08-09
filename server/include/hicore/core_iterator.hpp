#ifndef CORE_ITERATOR_HPP_
#define CORE_ITERATOR_HPP_

#include "core.hpp"
#include "mem.hpp"
#include <sstream>
#include <string>

namespace HiGIS::Core {

int64_t NextObject(size_t start, uint64_t *content, uint64_t size);

/** \brief Struct to represent a Polygon, use with PolygonAccessor */
typedef struct PolygonPtr {
    uint64_t *begin;
    uint64_t part_count;
    uint64_t *parts;
    uint64_t *rings;
    uint64_t *data;
} PolygonPtr;

class IWalker {
 public:
    virtual ~IWalker() {}
    virtual int64_t NextFeature() = 0;
    virtual int64_t NextFeature(uint64_t start) = 0;
};

/** \brief Mixin to add buffer-loop capabilities for Accessors, store current position in ptr */
template <typename AccessorT, typename PtrT> class Walker : public AccessorT, public IWalker {
 private:
    uint64_t *buffer_;
    uint64_t size_;
    void NextFeatureDirect(uint64_t index);

 public:
    Walker(uint64_t *buffer, uint64_t size, uint64_t start)
        : AccessorT(PtrT()), buffer_(buffer), size_(size) {
        NextFeature(start);
    }
    Walker(uint64_t *buffer, uint64_t size) : Walker(buffer, size, 0) {}
    Walker(GeoBuffer &buffer, uint64_t start) : Walker(buffer.Data(), buffer.Size(), start) {}
    Walker(GeoBuffer &buffer) : Walker(buffer, 0) {}
    /** \brief Next feature according to current ptr */
    int64_t NextFeature() {
        uint64_t index = AccessorT::ptr_.begin + AccessorT::Length() - buffer_;
        if (index < size_) {
            NextFeatureDirect(index);
            return index;
        } else {
            return -1;
        }
    }
    /**
     * \brief Next feature given start position
     *
     * @param start start position, as the index of buffer
     */
    int64_t NextFeature(uint64_t start) {
        auto index = NextObject(start, buffer_, size_);
        if (index >= 0) {
            NextFeatureDirect(index);
        }
        return index;
    }
};

/** \brief Class to access content of a polygon, like parts, rings, etc. */
class PolygonAccessor {
 protected:
    PolygonPtr ptr_;

 public:
    PolygonAccessor(const PolygonPtr &ptr) : ptr_(ptr) {}
    // light weight entity representations:
    // 1. polygon by its index
    // 2. part by its global index in l2
    size_t PartCount() { return ptr_.part_count; }
    size_t Part(size_t i) { return i; } // syntax sugar to be consistent
    /** \brief Ring count of the polygon part */
    size_t RingCount(size_t part) { return ptr_.parts[part + 1] - ptr_.parts[part]; }
    /** \brief The i-th ring of the part */
    size_t Ring(size_t part, size_t i) { return ptr_.parts[part] + i; }
    /** \brief The i-th ring of the whole, for quickly loop all rings */
    size_t Ring(size_t i) { return i; }
    /** \brief The node count the ring */
    size_t NodeCount(size_t ring) { return ptr_.rings[ring + 1] - ptr_.rings[ring]; }
    /** \brief Ring count in total */
    size_t RingCount() { return ptr_.parts[ptr_.part_count]; } // all ring count
    /** \brief Node count in total */
    size_t NodeCount() { return ptr_.rings[RingCount()]; } // all node count
    /** \brief All nodes */
    double *Nodes() { return reinterpret_cast<double *>(ptr_.data); }
    /** \brief Nodes of the ring */
    double *Nodes(size_t ring) {
        return reinterpret_cast<double *>(ptr_.data + ptr_.rings[ring] * 2);
    }
    double *Nodes(size_t ring, size_t start) {
        return reinterpret_cast<double *>(ptr_.data + ptr_.rings[ring] * 2) + start * 2;
    }
    /** \brief X of the i-th node of nodes */
    double X(double *nodes, size_t i) { return nodes[i * 2]; }
    /** \brief Y of the i-th node of nodes */
    double Y(double *nodes, size_t i) { return nodes[i * 2 + 1]; }
    PolygonPtr &Ptr() { return ptr_; }
    uint64_t Id() { return *ptr_.begin; }
    // (id, polygon count) = 2
    /** \brief The underlying memory consumed by the data, unit is 8 bytes */
    size_t Length();
};

/** \brief Class to loop a buffer and access content of a polygon */
using PolygonIterator = Walker<PolygonAccessor, PolygonPtr>;

/** \brief Struct to represent a LineString, use with LineStringAccessor */
typedef struct LineStringPtr {
    uint64_t *begin;
    uint64_t line_count;
    uint64_t *lines;
    uint64_t *data;
} LineStringPtr;

/** \brief Class to access content of a linestring, like lines, nodes, etc. */
class LineStringAccessor {
 protected:
    LineStringPtr ptr_;

 public:
    LineStringAccessor(const LineStringPtr &ptr) : ptr_(ptr) {}
    // light weight line representations by its index
    /** \brief Line count for MultiLineString */
    size_t LineCount() { return ptr_.line_count; }
    /** \brief The i-th line */
    size_t Line(size_t i) { return i; } // syntax sugar to be consistent
    /** \brief Node count of the line */
    size_t NodeCount(size_t line) { return ptr_.lines[line + 1] - ptr_.lines[line]; }
    /** \brief Node count in total */
    size_t NodeCount() { return ptr_.lines[LineCount()]; } // all node count
    /** \brief All nodes */
    double *Nodes() { return reinterpret_cast<double *>(ptr_.data); }
    /** \brief Nodes of the line */
    double *Nodes(size_t line) {
        return reinterpret_cast<double *>(ptr_.data + ptr_.lines[line] * 2);
    }
    double *Nodes(size_t ring, size_t start) {
        return reinterpret_cast<double *>(ptr_.data + ptr_.lines[ring] * 2) + start * 2;
    }
    double X(double *nodes, size_t i) { return nodes[i * 2]; }
    double Y(double *nodes, size_t i) { return nodes[i * 2 + 1]; }
    LineStringPtr &Ptr() { return ptr_; }
    uint64_t Id() { return *ptr_.begin; }
    // (id, line count) = 2
    /** \brief The underlying memory consumed by the data, unit is 8 bytes */
    size_t Length();
};

using LineStringIterator = Walker<LineStringAccessor, LineStringPtr>;

/** \brief struct to represent a MultiPoint, use with MultiPointAccessor */
struct MultiPointPtr {
    uint64_t *begin;
    uint64_t point_count;
    uint64_t *data;
};

/** \brief class to access content of a multipoint */
class MultiPointAccessor {
 protected:
    MultiPointPtr ptr_;

 public:
    MultiPointAccessor(const MultiPointPtr &ptr) : ptr_(ptr) {}
    size_t PointCount() { return ptr_.point_count; }
    double *Points() { return reinterpret_cast<double *>(ptr_.data); }
    double X(size_t i) { return Points()[i * 2]; }
    double Y(size_t i) { return Points()[i * 2 + 1]; }
    MultiPointPtr &Ptr() { return ptr_; }
    uint64_t Id() { return *ptr_.begin; }
    /** \brief The underlying memory consumed by the data, unit is 8 bytes */
    size_t Length() { return 2 + ptr_.point_count * 2 + 1; } // (id, point count) = 2
};

using MultiPointIterator = Walker<MultiPointAccessor, MultiPointPtr>;

/** \brief struct to represent a Point, use with PointAccessor */
typedef struct PointPtr {
    uint64_t *begin;
    uint64_t *data;
} PointPtr;

/** \brief class to access a point */
class PointAccessor {
 protected:
    PointPtr ptr_;

 public:
    PointAccessor(const PointPtr &ptr) : ptr_(ptr) {}
    double X() { return Coords()[0]; }
    double Y() { return Coords()[1]; }
    double *Coords() { return reinterpret_cast<double *>(ptr_.data); }
    bool Empty() { return *ptr_.data == TERMINATOR; }
    PointPtr &Ptr() { return ptr_; }
    uint64_t Id() { return *ptr_.begin; }
    /** \brief The underlying memory consumed by the data, (id, x, y, terminator) = 4 */
    size_t Length() { return Empty() ? 2 : 4; }
};
using PointIterator = Walker<PointAccessor, PointPtr>;

std::unique_ptr<IWalker> CreateWalker(GeoData &data);

// for human-readable printing, also illustrate usage

/** \ingroup core_iterator */
std::string Wkt(PointAccessor &iter);
/** \ingroup core_iterator */
std::string Wkt(MultiPointAccessor &iter);
/** \ingroup core_iterator */
std::string Wkt(PolygonAccessor &iter);
/** \ingroup core_iterator */
std::string Wkt(LineStringAccessor &iter);

} // namespace HiGIS::Core

#endif // CORE_ITERATOR_HPP_

/**
 * \defgroup core_iterator Core Iterator
 * \brief provides classes to iterator and access features' geometries, including Ptrs,
 * Accessors and Iterators.
 *
 *
 * *Iterators*
 *
 * The main usage pattern for this module is to use Iterators to loop a feature pool and access the
 * detail of each feature. Iterators mix Accessor APIs and Walker APIs.
 *
 * Build a Wkt function for single features using Accessor APIs:
 *
 * \snippet core_iterator.cpp PolygonComponents
 *
 * \code
 * PolygonIterator iter(buffer);  // buffer: GeoBuffer
 * do {
 *     std::cout << Wkt(iter) << std::endl;
 * } while (iter.NextFeature() >= 0);
 * \endcode
 *
 * \see HiGIS::Core::Walker, HiGIS::Core::PolygonAccessor
 *
 * Sometimes we may deal with irregular data, in that case we may want to manually construct *Ptrs*
 * and use *Accessors* to access it.
 *
 * Ptrs hold data necessary to describe a feature. Accessors access the memory structure of a
 * feature, basically it's only a thin wrapper providing syntax sugar on Ptrs. The polygon, ring,
 * line, etc. returned by functions like `Ring(size_t)`, `Line(size_t)` are all indices by which the
 * subsequent method can find the details of the returned entities. This design prevents small
 * object allocations and enables some quite flexible usage patterns. Accessors have the same API as
 * Iterators except for the constructors and NextFeature methods (In fact, Iterators mix-in
 * Accessors).
 */

/**
 * \class HiGIS::Core::Walker
 * \ingroup core_iterator
 *
 * Usage example: \code
 * // return -1 if reach the end
 * walker.NextFeature(); // next feature according to current ptr
 * walker.NextFeature(200); // next feature searching from the-200-th buffer unit
 * \endcode
 *
 * \see core_iterator
 */

/**
 * \class HiGIS::Core::PolygonAccessor
 * \ingroup core_iterator
 * \snippet core_iterator.cpp PolygonComponents
 *
 * \class HiGIS::Core::LineStringAccessor
 * \ingroup core_iterator
 * \snippet core_iterator.cpp LineStringComponents
 * 
 * \class HiGIS::Core::MultiPointAccessor
 * \ingroup core_iterator
 * \snippet core_iterator.cpp MultiPointComponents
 * 
 * \class HiGIS::Core::PointAccessor
 * \ingroup core_iterator
 * \snippet core_iterator.cpp PointComponents
 */

/**
 * \class HiGIS::Core::PolygonIterator
 * \ingroup core_iterator
 * The class can use both APIs from PolygonAccessor and Walker
 * 
 * \copydoc PolygonAccessor
 * \code 
 * do {
 *     // ...
 * } while(iter.NextFeature() >= 0)
 * \endcode
 */

/** 
 * \class HiGIS::Core::LineStringIterator
 * \ingroup core_iterator
 * The class can use both APIs from LineStringAccessor and Walker
 * 
 * \copydoc LineStringAccessor
 * \code 
 * do {
 *     // ...
 * } while(iter.NextFeature() >= 0)
 * \endcode
 */

/** \class HiGIS::Core::MultiPointIterator
 * \ingroup core_iterator
 * The class can use both APIs from MultiPointAccessor and Walker
 * 
 * Accessor:
 * \copydoc MultiPointAccessor
 * \code 
 * do {
 *     // ...
 * } while(iter.NextFeature() >= 0)
 * \endcode
 */

/**
 * \class HiGIS::Core::PointIterator
 * \ingroup core_iterator
 * The class can use both APIs from PointAccessor and Walker
 * 
 * Accessor:
 * \copydoc PointAccessor
 * \code 
 * do {
 *     // ...
 * } while(iter.NextFeature() >= 0)
 * \endcode
 */