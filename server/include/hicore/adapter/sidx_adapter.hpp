#ifndef LIBSPATIALINDEX_ADAPTER_HPP
#define LIBSPATIALINDEX_ADAPTER_HPP

#include <spatialindex/SpatialIndex.h>

#include "core_iterator.hpp"
#include "geo.hpp"
#include "mem.hpp"
#include "util.hpp"

namespace HiGIS::Adapter {
using Core::Extent;
using SpatialIndex::Point;
using SpatialIndex::Region;

Extent RegionToExtent(const Region *r) {
    return Extent{r->m_pLow[0], r->m_pLow[1], r->m_pHigh[0], r->m_pHigh[1]};
}

class SIPoint : public virtual SpatialIndex::IShape {
 public:
    SIPoint(Core::PointIterator &iter) : iter_(iter) {}
    ~SIPoint() = default;
    //
    // ISerializable interface
    //
    uint32_t getByteArraySize() override;
    void loadFromByteArray(const uint8_t *data) override;
    void storeToByteArray(uint8_t **data, uint32_t &length) override;

    //
    // IShape interface
    //
    bool intersectsShape(const IShape &in) const override;
    bool containsShape(const IShape &in) const override;
    bool touchesShape(const IShape &in) const override;
    void getCenter(Point &out) const override;
    uint32_t getDimension() const override;
    void getMBR(Region &out) const override;
    double getArea() const override;
    double getMinimumDistance(const IShape &in) const override;

 private:
    Core::PointIterator &iter_;
};

inline bool SIPoint::intersectsShape(const SpatialIndex::IShape &in) const {
    const Region *pr = dynamic_cast<const Region *>(&in);
    if (pr != nullptr) {
        return ContainsPoint(RegionToExtent(pr), iter_.X(), iter_.Y());
    }
    throw Tools::IllegalStateException("Point::intersectsShape: Not implemented yet!");
}
inline bool SIPoint::containsShape(const SpatialIndex::IShape &) const { return false; }
inline bool SIPoint::touchesShape(const SpatialIndex::IShape &in) const {
    const Region *pr = dynamic_cast<const Region *>(&in);
    if (pr != nullptr) {
        return TouchesPoint(RegionToExtent(pr), iter_.X(), iter_.Y());
    }
    throw Tools::IllegalStateException("Point::touchesShape: Not implemented yet!");
}
inline void SIPoint::getCenter(Point &out) const { out = Point(iter_.Coords(), 2); }
inline uint32_t SIPoint::getDimension() const { return 2; }
inline void SIPoint::getMBR(Region &out) const { out = Region(iter_.Coords(), iter_.Coords(), 2); }
inline double SIPoint::getArea() const { return 0; }
inline double SIPoint::getMinimumDistance(const SpatialIndex::IShape &in) const {
    auto ppt = dynamic_cast<const SIPoint *>(&in);
    if (ppt != nullptr) {
        return Distance(iter_.X(), iter_.Y(), ppt->iter_.X(), ppt->iter_.Y());
    }
    const Region *pr = dynamic_cast<const Region *>(&in);
    if (pr != nullptr) {
        return GetMinDistPoint(RegionToExtent(pr), iter_.X(), iter_.Y());
    }
    throw Tools::IllegalStateException("Point::getMinimumDistance: Not implemented yet!");
}
inline uint32_t SIPoint::getByteArraySize() { return iter_.Length() * 8; }
inline void SIPoint::loadFromByteArray(const uint8_t *data) {
    iter_ = Core::PointIterator(const_cast<uint64_t *>(reinterpret_cast<const uint64_t *>(data)),
                                getByteArraySize() / 8);
}
inline void SIPoint::storeToByteArray(uint8_t **data, uint32_t &length) {
    *data = reinterpret_cast<uint8_t *>(iter_.Ptr().begin);
    length = getByteArraySize();
}

} // namespace HiGIS::Adapter

#endif