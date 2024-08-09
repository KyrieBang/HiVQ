#ifndef BINARY_HPP_
#define BINARY_HPP_

#include "mem.hpp"
#include <memory>
#include <stdlib.h>
#include <string>

namespace HiGIS {

template <typename T, int n> struct chunk_t { T data[n]; };

using endian_t = bool;
const endian_t BE = true;
const endian_t LE = false;

endian_t SysEndian(void);

void SwapEndian(uint64_t &value);
void SwapEndianBatch(uint8_t *array, size_t bunk_size, size_t num);
int GetInt4(endian_t endian, const uint8_t *cursor);
double GetDouble8(bool same_endian, const uint8_t *cursor);
int HexToBytes(uint8_t *out, const std::string &hex);
} // namespace HiGIS

#endif // BINARY_HPP_