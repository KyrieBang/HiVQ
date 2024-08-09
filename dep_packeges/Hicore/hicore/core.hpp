#ifndef CORE_HPP_
#define CORE_HPP_

#include "mem.hpp"
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

namespace HiGIS::Core {

/**
 * \brief Mix multiple structs by multi-inherit
 *
 * Mix provides a way to combine data of several structs without the need of verbose redirections.
 * Mix<A, B, C> will inherit A, B, C and thus have all members of A, B, C.
 *
 * Polymorphism is somehow permitted. Rules to pass Mix as arguments:
 *
 * 1. Allow `Mix<A, B, C> -> Mix<A, B>`
 * 2. Allow `Mix<A, B, C> -> &A or Mix<A, B, C> -> A`
 * 3. Do not allow `Mix<A, B, C> -> Mix<A, B>&`, use `Mix<A, B>(A&, B&)` to reconstruct
 */
template <typename... Is> struct Mix : public virtual Is... {
 public:
    Mix() : Is()... {}
    /**
     * \brief Allow construct from Mix with more Mixed structs, like from Mix<A, B, C> to Mix<A, B>.
     *
     * This enables implicit conversion when passing arguments (Polymorphism rule 1)
     */
    template <typename T> Mix(T &t) : Is(t)... {}
    /**
     * \brief Allow construct from Mixed components, like from A, B, C to Mix<A, B, C>.
     *
     * This enables flexible reconstruction of Mix (Polymorphism rule 3)
     */
    Mix(Is &... ts) : Is(ts)... {}
};

/**
 * \brief Terminator at the end of each feature
 * \ingroup core
 */
const uint64_t TERMINATOR = ~0LL;

enum GeoType {
    Unknown = -1,
    /** Also store simple polygon as a one-part multipolygon */
    MultiPolygon = 6,
    /** Also store simple linestring as a one-part multilinestring */
    MultiLineString = 5,
    MultiPoint = 4,
    Point = 1
};

using GeoBuffer = MemPool<uint64_t>;

struct Extent {
    double xmin = INFINITY;
    double ymin = INFINITY;
    double xmax = -INFINITY;
    double ymax = -INFINITY;
    double *Values() { return reinterpret_cast<double *>(this); }
};

// Stats
struct FeatureExtents {
    MemPool<Extent> extents;
};
using GeoStats = Mix<FeatureExtents>;

struct GeoSummary {
    uint64_t feature_count;
    GeoType feature_type;
    struct Extent extent;

    std::string srs;
    // GeoSummary() : feature_count(0), feature_type(GeoType::Unknown) {}
};

struct Index {
    MemPool<uint64_t> id_index;
};

template <typename TIndex> using GeoDataT = Mix<GeoBuffer, GeoSummary, GeoStats, TIndex>;
using GeoData = Mix<GeoBuffer, GeoSummary, GeoStats>;

} // namespace HiGIS::Core

#endif // CORE_HPP_

/**
 * \defgroup core Core
 * \brief Basic structures for HiGIS Core
 *
 * Generally, a geographical dataset consists of 4 possible parts:
 *
 * - Buffer for binary data
 * - Summary for overall description of the dataset, like feature count, total area, etc.
 * - Stats for descriptions of each feature, like extents, polygon areas, line lengths, etc.
 * - Index for spatial index for the dataset (not implemented yet)
 *
 * The module use Mix to provide a standard version of dataset named GeoData, consisting with a
 * standard version of Summary(GeoSummary), Stats(GeoStats)
 */

/**
 * \class HiGIS::Core::Mix
 * \ingroup core
 * \snippet test_core.cpp MixTest
 */

/**
 * \class GeoSummary
 * \ingroup core
 * \brief Basic definition for Summary
 */

/**
 * \class GeoBuffer
 * \ingroup core
 * \brief Buffer to store feature geometries
 *
 * The unit is universally 8 bytes instead of keeping variant units like 1, 4, 8 bytes, which is
 * space-consuming but is much easier to deal with Especially, it makes alignment possible so we can
 * scan from arbitrary points. Point coordinates as double are simply reinterprated and put into the
 * buffer, with its memory representation intact.
 */

/**
 * \class GeoStats
 * \ingroup core
 * \brief Basic definition for Stats
 */

/**
 * \class FeatureExtents
 * \ingroup core
 * \brief FeatureExtents Mixin for Stats
 */

/**
 * \class GeoData
 * \ingroup core
 * \brief Useful data definition for most cases
 */

/**
 * \enum GeoType
 * \ingroup core
 *
 * The types are basically from Shapefile. Polygon and MultiPolygon, LineString and MultiLineString
 * are not distinguished. The reason is that restricting polygons to have only one exterior ring is
 * not very useful. The non-multi versions are simply regarded as special cases.
 */