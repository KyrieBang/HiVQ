#ifndef INDEX_HPP
#define INDEX_HPP

#include <cstdint>
#include <sstream>
#include <string>

std::string CodeToStr(uint64_t code, int length) {
    uint64_t mask = 3;
    std::stringstream ss;
    for (size_t offset = 0; offset < length; offset += 2) {
        ss << ((code >> offset) & mask);
    }
    return ss.str();
}

namespace HiGIS {
class QuadTree {
 private:
    uint64_t node_max_;

 public:
    uint64_t count;
    uint64_t code;

    QuadTree(uint64_t node_max);
    ~QuadTree();
};

} // namespace HiGIS

#endif