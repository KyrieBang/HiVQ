#include "hicore/core.hpp"
#include "hicore/core_creator.hpp"
#include "hicore/core_iterator.hpp"
#include "hicore/shp.hpp"
#include "treeNode.hpp"

# ifndef BUILDINDEX_HPP_
# define BUILDINDEX_HPP_

using namespace HiGIS::IO;
using namespace HiGIS::Core;


void coordTran(double &x, double &y);
pair<pointNode*, short> pointIndex_shp(string data_path, string data_srs, FieldSet part_fields);
pair<linestringNode*, short> linestringIndex_shp(string data_path, string data_srs, FieldSet part_fields);
pair<polygonNode*, short> polygonIndex_shp(string data_path, string data_srs, FieldSet part_fields);

# endif