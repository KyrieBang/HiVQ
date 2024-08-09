#include "core.hpp"
#include "core_iterator.hpp"
#include <memory>
#include <string>

#ifndef SHP_HPP_
#define SHP_HPP_

namespace HiGIS::IO {

using Core::GeoBuffer;
using Core::GeoData;
using Core::GeoDataT;
using Core::GeoSummary;
using Core::GeoType;

/** \brief Disk loader, deciding how to map binary file into the memory and release */
class IDiskLoader {
 public:
    /**
     * \brief Load file into memory
     * \param path file path
     * \return pointer to the memory
     */
    virtual uint8_t *Load(const std::string &path) = 0;
    /** \brief release loaded data */
    virtual void Unload() = 0;
    virtual ~IDiskLoader() {}
    /** \brief the byte count of loaded data in memory */
    virtual uint64_t Size() { return size_; }

 protected:
    uint64_t size_;
};

/**
 * \brief Disk loader using mmap, slower than fstream loader in the first run, faster in the
 * second run
 * \ingroup shape
 */
class MmapLoader : public IDiskLoader {
 private:
    uint8_t *buffer_;
    int fd_;

 public:
    uint8_t *Load(const std::string &path);
    void Unload();
};

/**
 * \brief Disk loader using fstream, the default loader
 * \ingroup shape
 */
class FreadLoader : public IDiskLoader {
 private:
    uint8_t *buffer_;

 public:
    uint8_t *Load(const std::string &path);
    void Unload();
};

/**
 * \brief Disk loader using fstream, the default loader
 * \ingroup shape
 */
class ReadLoader : public IDiskLoader {
 private:
    uint8_t *buffer_;

 public:
    uint8_t *Load(const std::string &path);
    void Unload();
};

/**
 * \brief Disk loader using fstream, the default loader
 * \ingroup shape
 */
class FstreamLoader : public IDiskLoader {
 private:
    uint8_t *buffer_;

 public:
    uint8_t *Load(const std::string &path);
    void Unload();
};

/** \brief The meta section of a shapefile */
struct ShapeMeta {
    int32_t file_code;
    int32_t file_length;
    int32_t version;
    int32_t shape_type;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    // real file size, file_length may be misleading for very large files
    uint64_t size;
};

/** \brief binary of the components of a shapefile */
struct OutBuffers {
    MemPool<uint8_t> shp_out;
    MemPool<uint8_t> shx_out;
};

/** \brief Reader and writer of a shapefile */
class Shapefile {
 public:
    /** \brief Shapefile options to change the behaviors when read & write */
    struct Options {
        /** \brief The method to load file into memory */
        std::shared_ptr<IDiskLoader> disk_loader;
        /** \brief Whether to totally skip empty features */
        bool skip_empty;
        Options() {
            disk_loader = std::make_shared<ReadLoader>();
            skip_empty = false;
        }
    };

    Shapefile(Options opt)
        : options_(opt), sys_endian_(SysEndian()), sys_le_(sys_endian_ == LE),
          sys_be_(sys_endian_ == BE), temp_(std::make_shared<GeoBuffer>(2097152)) {}

    /** \brief Shapefile constructor with default options */
    Shapefile() : Shapefile(Options()) {}

    /** \brief Read a shapefile from path, equal to Load then Read memory */
    GeoData Read(const std::string &path);
    /**
     * \brief Read in-memory shapefile
     * \param buffer pointer to the memory to read
     * \param size The size of the memory in bytes
     * \return The result data
     */
    GeoData Read(uint8_t *buffer, uint64_t size = 0);
    /**
     * \brief Read the header from in-memory shapefile, the result is a precise clone of the
     * shapefile header
     *
     * \param buffer pointer to the memory to read
     * \param file_size The size of the * memory in bytes
     * \return The meta data
     */
    ShapeMeta ReadShpMeta(uint8_t *buffer, uint64_t file_size = 0);
    /**
     * \brief Read the meta data from in-memory shapefile, the result is the
     * organized meta data
     *
     * \param out The summary object to fill
     * \param buffer pointer to the memory to read
     * \return The cursor after reading meta
     */
    uint8_t *ReadMeta(GeoSummary &out, uint8_t *buffer);
    /**
     * \brief Read the next feature in memory
     *
     * \param out The buffer object to fill
     * \param id The id of the result feature, passed in so that it can be different from the
     * original id
     * \param cursor The pointer to the current location of memory
     * \param type The geometry type
     * \return The cursor after reading the feature
     */
    uint8_t *ReadFeature(GeoBuffer &out, uint64_t id, uint8_t *cursor, Core::GeoType type);
    /**
     * \brief Read the extent in memory
     *
     * \param out The extent output buffer to append
     * \param cursor The pointer to the current location of memory
     * \return The cursor after reading the extent
     */
    void ReadExtent(MemPool<Core::Extent> &out, uint8_t *cursor);
    /** \brief Update the layer summary */
    void UpdateSummary(GeoSummary &out, uint8_t *cursor);
    /**
     * \brief Write data to a file
     *
     * \param data The data to write
     * \param path The output file path
     */
    void Write(GeoData &data, const std::string &path);
    /** \brief Write to memory */
    OutBuffers Write(GeoData &data);
    /**
     * \brief The path of auxiliary files of a shapefile, like .shx, .dbf, .prf etc.
     *
     * \param shp_path Path of the shp file
     * \param ext The file extension
     * \return The auxiliary file path
     */
    std::string AuxPath(std::string shp_path, const std::string &ext);

 private:
    Options options_;
    endian_t sys_endian_;
    endian_t sys_le_;
    endian_t sys_be_;
    std::shared_ptr<GeoBuffer> temp_;

    template <typename IterT>
    void WriteFeatures(MemPool<uint8_t> &shp_out, MemPool<uint8_t> &shx_out, GeoData &data);
    template <typename IterT> void WriteFeature(MemPool<uint8_t> &out, IterT &iter, GeoData &data);
    uint8_t *ReadPolygon(GeoBuffer &out, uint8_t *buffer);
    uint8_t *ReadLineString(GeoBuffer &out, uint8_t *buffer);
    uint8_t *ReadMultiPoint(GeoBuffer &out, uint8_t *buffer);
    uint8_t *ReadPoint(GeoBuffer &out, uint8_t *buffer);
    void WriteMeta(ShapeMeta &, MemPool<uint8_t> &buffer);
    std::string ReadSrs(std::string base_path);
    void WriteSrs(std::string base_path, const std::string &srs);
};

// the function is inline here deliberately for reference
/**
 * \brief Read in-memory shapefile, also serves as a reference to write read functions for custom
 * data structure
 *
 * \param shp The shapefile reader
 * \param buffer Pointer to the memory to read
 * \param file_size The size of the memory in bytes
 * \return The result data
 */
inline GeoData Read(Shapefile &shp, uint8_t *buffer, uint64_t file_size = 0) {
    GeoData data;
    auto shp_meta = shp.ReadShpMeta(buffer, file_size);
    data.Grow(shp_meta.size / 6);         // heuristic, about 1.3 * original file size
    data.extents.Grow(shp_meta.size / 8); // heuristic, 1 * original file size

    auto cursor = shp.ReadMeta(data, buffer);
    while (cursor < buffer + shp_meta.size) {
        shp.UpdateSummary(data, cursor);
        shp.ReadExtent(data.extents, cursor);
        cursor = shp.ReadFeature(data, data.feature_count, cursor, data.feature_type);
    }

    return data;
}

inline GeoDataT<Core::Index> ReadWithIdIndex(Shapefile &shp, uint8_t *buffer,
                                             uint64_t file_size = 0) {
    auto data = Read(shp, buffer, file_size);
    Core::Index index = {MemPool<uint64_t>(data.feature_count)};
    auto walker = Core::CreateWalker(data);
    int64_t i = 0;
    do {
        index.id_index.Push(i);
        i = walker->NextFeature();
    } while (i >= 0);
    return GeoDataT<Core::Index>(data, data, data, index);
}

} // namespace HiGIS::IO

#endif // SHP_HPP_

/**
 * \defgroup shape Shapefile
 * \brief Read and write shapefiles
 *
 * Main class: HiGIS::IO::Shapefile
 * \copydoc HiGIS::IO::Shapefile
 *
 * \class HiGIS::IO::Shapefile
 * \ingroup shape
 *
 * Usage example:
 * \code
 *     // Read
 *     auto out = Shapefile().Read("~/data/test/polygon.shp");
 *     // Write
 *     Shapefile().Write(out, "~/data/test/polygon_copy.shp");
 *
 *     // Use mmap to read
 *     Shapefile::Options opts = { std::make_shared<MmapLoader>() };
 *     auto out = Shapefile(opts).Read("~/data/test/polygon.shp");
 * \endcode

 */