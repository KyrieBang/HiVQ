#ifndef GEO_HPP_
#define GEO_HPP_
#include "core.hpp"

namespace HiGIS {

void UpdateExtent(Core::Extent &ext, const Core::Extent &other);
double AreaSign(int n, double xy[]);
// clockwise
inline bool IsOuterRing(int n, double xy[]) { return AreaSign(n, xy) < 0; }
inline double Area(int n, double xy[]) { return AreaSign(n, xy) / 2.0; }
// inline double Area(const Core::Extent &extent);
bool PointInRing(double x, double y, double *xy, int n);
bool ContainsPoint(const Core::Extent &extent, double x, double y);
// bool ContainsRegion(const Core::Extent &extent, const Core::Extent &other);
bool TouchesPoint(const Core::Extent &extent, double x, double y);
// bool TouchesRegion(const Core::Extent &extent, double x, double y);
// bool IntersectsSegment(const Core::Extent &extent, double x1, double y1, double x2, double y2);
// bool IntersectsRegion(const Core::Extent &extent, const Core::Extent &other);
// bool GetCenter(const Core::Extent &extent);
double GetMinDistPoint(const Core::Extent &extent, double x, double y);
// double GetMinDistExtent(const Core::Extent &extent, const Core::Extent &other);
double Distance(double x1, double y1, double x2, double y2);

} // namespace HiGIS // namespace HiGIS

#endif
