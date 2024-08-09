#ifndef MEM_IMPL_HPP_
#define MEM_IMPL_HPP_

#include "binary.hpp"
#include "mem.hpp"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>

namespace HiGIS {

template <typename T>
MemPool<T>::MemPool(uint64_t chunk_size) : chunk_size_(chunk_size), n_(0), capacity_(0) {}

template <typename T> MemPool<T>::MemPool() : MemPool(MEMPOOL_CHUNK_SIZE) {}
template <typename T> MemPool<T>::~MemPool() {}

template <typename T> void MemPool<T>::Grow(uint64_t size) {
    auto new_size = (capacity_ + size) * sizeof(T);
    auto new_buffer_ = std::shared_ptr<T>((T *)malloc(new_size), free);
    if (new_buffer_ == nullptr) {
        std::cerr << "Failed to allocate size " + std::to_string(new_size) << std::endl;
        throw std::bad_alloc();
    }
    memcpy(new_buffer_.get(), buffer_.get(), capacity_ * sizeof(T));
    buffer_ = new_buffer_;
    capacity_ += size;
}

template <typename T> void MemPool<T>::Grow() { Grow(chunk_size_); }

template <typename T> void MemPool<T>::Push(T item) {
    if (n_ >= capacity_) {
        Grow();
    }
    Data()[n_] = item;
    n_++;
}

template <typename T> void MemPool<T>::Copy(T *items, size_t n) {
    Reserve(n);
    memcpy(Data() + n_ - n, items, n * sizeof(T));
}

template <typename T> void MemPool<T>::Reserve(uint64_t n) {
    uint64_t new_n = n_ + n;
    if (new_n > capacity_) {
        Grow(((new_n - capacity_) / chunk_size_ + 1) * chunk_size_);
    }
    n_ = new_n;
}

template <typename T> T *MemPool<T>::Data() { return buffer_.get(); }
template <typename T> T *MemPool<T>::Cursor() { return buffer_.get() + n_; }
template <typename T> void MemPool<T>::Clear() { n_ = 0; }
template <typename T> void MemPool<T>::Clear(uint64_t n) { n_ = n; }

template <typename T> void WriteBytes(bool same_endian, T value, MemPool<uint8_t> &out) {
    auto n = sizeof(T);
    auto buf = reinterpret_cast<uint8_t *>(&value);
    for (size_t i = 0; i < n; i++) {
        if (same_endian) {
            out.Push(buf[i]);
        } else {
            out.Push(buf[n - 1 - i]);
        }
    }
}

template <typename T>
void WriteBytes(bool same_endian, T value, MemPool<uint8_t> &out, uint64_t index) {
    auto n = sizeof(T);
    auto buf = reinterpret_cast<uint8_t *>(&value);
    for (size_t i = 0; i < n; i++) {
        if (same_endian) {
            out.Data()[index + i] = buf[i];
        } else {
            out.Data()[index + i] = buf[n - 1 - i];
        }
    }
}

template <typename T> std::string MemPool<T>::Hex() {
    auto buffer = reinterpret_cast<uint8_t *>(Data());
    std::stringstream ss;
    for (size_t i = 0; i < n_ * sizeof(T); i++) {
        ss << std::hex << std::setfill('0') << std::setw(2) << std::uppercase
           << (int)buffer[i]; // has to be int since uint will print char anyway
    }

    return ss.str();
}
template <typename T> T MemPool<T>::operator[](int64_t i) const {
    if (i < 0)
        return Data()[Size() + i];
    return Data()[i];
}
template <typename T> T &MemPool<T>::operator[](int64_t i) {
    if (i < 0)
        return Data()[Size() + i];
    return Data()[i];
}

} // namespace HiGIS

#endif // MEM_IMPL_HPP_