#include "core.hpp"
#include <mpi.h>
#include <stddef.h>

namespace HiGIS {

using Core::GeoSummary;

// T should be child of GeoSummary
template <typename T> class MPIHelper {
 public:
    MPI_Datatype GEOSUMMARY_TYPE;

    MPIHelper() {
        static_assert(std::is_base_of<GeoSummary, T>::value, "T must inherit from GeoSummary");
        const int nfield = 6;
        const int nitems = nfield;
        int blocklengths[nfield] = {1, 1, 1, 1, 1, 1};
        MPI_Datatype types[nfield] = {
            MPI_UNSIGNED_LONG_LONG, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
        MPI_Aint offsets[nfield] = {offsetof(T, feature_count), offsetof(T, feature_type),
                                    offsetof(T, extent.xmin),   offsetof(T, extent.ymin),
                                    offsetof(T, extent.xmax),   offsetof(T, extent.ymax)};

        MPI_Type_create_struct(nitems, blocklengths, offsets, types, &GEOSUMMARY_TYPE);
        MPI_Type_commit(&GEOSUMMARY_TYPE);
    }

    ~MPIHelper() { MPI_Type_free(&GEOSUMMARY_TYPE); }
};

} // namespace HiGIS