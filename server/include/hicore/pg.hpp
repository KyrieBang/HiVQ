#include "core.hpp"
#include <libpq-fe.h>
#include <memory>
#include <string>

#ifndef PG_HPP_
#define PG_HPP_

namespace HiGIS::IO {

using std::string;

class PostGIS;

struct PGresultPtr {
    PGresult *data;
    PGconn *conn_;
    PGresultPtr(PostGIS &pg, const std::string &table);
    ~PGresultPtr();
};

class PostGIS {
 public:
    // password is not passed by options, use ~/.pgpass file
    struct Options {
        string host = "localhost";
        int port = 5432;
        string user = "postgres";
        string dbname;
        string geometry_column = "wkb_geometry";
    };

    PostGIS(Options opt);
    ~PostGIS();

    void ReadFeature(Core::GeoBuffer &out, PGresult *data, int type, uint64_t i);
    int ReadMeta(Core::GeoSummary &out, PGresult *data);

 private:
    Options options_;
    PGconn *conn_;
    friend struct PGresultPtr;
};

inline Core::GeoData Read(PostGIS &pg, const std::string &table) {
    Core::GeoData out;
    PGresultPtr res_ptr(pg, table);
    auto type = pg.ReadMeta(out, res_ptr.data);
    for (uint64_t i = 0; i < out.feature_count; i++) {
        pg.ReadFeature(out, res_ptr.data, type, i);
    }
    return out;
}

} // namespace HiGIS::IO

#endif // PG_HPP_
