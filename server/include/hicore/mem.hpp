#ifndef MEM_HPP_
#define MEM_HPP_

#include "binary.hpp"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdlib.h>
#include <string>

namespace HiGIS {

const uint64_t MEMPOOL_CHUNK_SIZE = 4194304; // for size 8, take 32MB

/**
 * \brief A memory pool with continous memory on heap, will grow if not large enough
 */
template <typename T> class MemPool {
 private:
    std::shared_ptr<T> buffer_;
    uint64_t chunk_size_;
    uint64_t n_;
    uint64_t capacity_;

 public:
    /** \brief = MemPool(MEMPOOL_CHUNK_SIZE) */
    MemPool();
    /**
     * \brief Specifying the space expand every time the pool is too small
     * @param chunk_size allocate current_size + chunk_size and migrate the data when grow
     */
    MemPool(uint64_t chunk_size);
    /** \brief Deconstructor, will free the allocated memory */
    ~MemPool();

    /** \brief push an item into the pool, grow if necessary */
    void Push(T item);
    /** \brief copy a memory block to the pool, grow if necessary */
    void Copy(T *items, size_t n);
    /** \brief reserve empty space of n * (item size), grow if necessary */
    void Reserve(uint64_t n);
    /** \brief clear the pool (will not shrink the size) */
    void Clear();
    /** \brief clear the pool after n (will not shrink the size) */
    void Clear(uint64_t n);
    /** \brief grow the pool to size + chunk_size, should no call manually in most case */
    void Grow();
    /** \brief grow the pool to size + size_to_grow, should no call manually in most case */
    void Grow(uint64_t size_to_grow);
    uint64_t Size() { return n_; }
    /** \brief row pointer to the underlying memory, return value may be invalid if grow happends
     * afterwards */
    T *Data();
    /** \brief raw pointer to the underlying memory at current position to push, return value may be
     * invalid if grow happens afterwards */
    T *Cursor();
    /** \brief Use pool[i] to get element, -i to get in reverse order */
    T operator[](int64_t i) const;
    /** \brief Use pool[i] = value to set element, -i to get in reverse order the i-th position
     * should already has data */
    T &operator[](int64_t i);
    /** \brief print hex presentation of the data in pool */
    std::string Hex();
};

} // namespace HiGIS
#include "mem_impl.hpp"
#endif // MEM_HPP_

/**
 * \class HiGIS::MemPool
 * \ingroup core
 *
 * MemPool is a convenient data structure for storing items in continuous memory. The space is only
 * allocated when necessary so operations will be much faster than allocating many small objects.
 *
 * MemPool acts like vector, but the allocating strategy is linear, which is more reasonable for
 * large data. The API style, unlike vector, is not heavily templated, but is often in a C-like
 * manner.
 *
 * MemPool is append-only.
 *
 * Usage example: \code
 * // Declare a new MemPool
 * HiGIS::MemPool<int> pool;
 * // Or a new MemPool providing a chunk_size of 256KB (default is 256MB)
 * HiGIS::MemPool<int> pool(256*1024);
 * \endcode
 *
 * Test cases to illustrate full APIs
 * \snippet test_mem.cpp MemPoolTest
 */
