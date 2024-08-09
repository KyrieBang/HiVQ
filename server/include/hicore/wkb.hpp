#ifndef WKB_HPP_
#define WKB_HPP_

#include "core.hpp"

#define POLYGON_TYPE 3
#define LINESTRING_TYPE 2

namespace HiGIS::IO {

Core::GeoType ReadWkb(Core::GeoBuffer &out, uint8_t *source, int type);
int ReadWkbType(uint8_t *source);
Core::GeoType GetGeoType(int code);

} // namespace HiGIS::IO

#endif // WKB_HPP_
