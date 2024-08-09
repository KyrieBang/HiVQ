#include <iostream>
#include <boost/geometry.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <sys/time.h>
#include <map>
#include <shapefil.h>

using namespace std;
namespace bg	= boost::geometry;
namespace bgi	= boost::geometry::index;
namespace bgm	= boost::geometry::model;

#ifndef _TREENODE_HPP_
#define _TREENODE_HPP_


enum nodeType{
    ROOT,
    LU,
    RU,
    LB,
    RB,
};


struct pointAttri
{
    size_t id;
    map<string, string> field_string;
};

struct linestringAttri
{
    size_t id;
    map<string, string> field_string;
};

struct polygonAttri
{
    size_t id;
    map<string, string> field_string;
};

struct FieldSet
{
    vector<string> fields_classify;
    vector<string> fields_continuous;
};



typedef bgm::d2::point_xy<double> POINT;
typedef bgm::box<POINT> Box;
typedef pair<POINT, pointAttri> POINT_ITEM;
typedef bgi::rtree<POINT_ITEM, bgi::quadratic<8>> pointRtree;
typedef bgm::segment<POINT>	SEG;
typedef vector<double> LINE;
typedef bgm::linestring<POINT> LINESTRING;
typedef tuple<Box, LINE, linestringAttri> LINESTRING_ITEM;
// typedef tuple<Box, string, linestringAttri> LINESTRING;
// typedef pair<Box, int> LINESTRING;
typedef bgi::rtree<LINESTRING_ITEM, bgi::quadratic<8>> linestringRtree;
typedef vector<double> RING;
typedef bgm::polygon<POINT> POLYGON;
typedef tuple<Box, RING, polygonAttri> POLYGON_ITEM;
typedef bgi::rtree<POLYGON_ITEM, bgi::quadratic<8>> polygonRtree;

// typedef map<string, pair<int, double>> COUNTER;
typedef map<string, map<string, int>> POINT_ATTRI_NODE;
typedef map<string, map<string, double>> LINE_ATTRI_NODE;
typedef map<string, map<string, double>> POLYGON_ATTRI_NODE;


class pointNode{
    public:
        // node type
        nodeType node_type;
        // node BBox
        double node_x; 
        double node_y;
        double node_w;
        // node attributes
        unsigned short node_level;
        unsigned int node_fcount;
        // float ANN;
        POINT_ATTRI_NODE attri_map;
        // node ptrs
        pointNode *parent;
        pointNode *LUNode;
        pointNode *RUNode;
        pointNode *LBNode;
        pointNode *RBNode;
        pointRtree Rtree;
        // node functions
        pointNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, pointNode *_parent);
        ~pointNode();
        void insertObj(double obj_x, double obj_y, short switch_level, pointAttri attri);
        void saveIndex(ofstream &outFile);
        pointNode *getMinWNode(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level);
        void countField_single(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, map<string, int> &count_map);
        void countField_classified(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, map<string, int> &count_map);
        void countField_ruled(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, vector<string> rules, map<string, int> &count_map);
        void getMaxFcount(int z, int &max_fcount);
        void getMaxFcount_range(int z, int &max_fcount, double box[]);
};

class linestringNode{
    public:
        // node type
        nodeType node_type;
        // node BBox
        double node_x; 
        double node_y;
        double node_w;
        // node attributes
        unsigned short node_level;
        double node_flen;
        // float ANN;
        LINE_ATTRI_NODE attri_map;
        // node ptrs
        linestringNode *parent;
        linestringNode *LUNode;
        linestringNode *RUNode;
        linestringNode *LBNode;
        linestringNode *RBNode;
        linestringRtree Rtree;
        // node functions
        linestringNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, linestringNode *_parent);
        ~linestringNode();
        linestringNode *insertTMBB(double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax, int seg_count, double len, short switch_level, linestringAttri attri);
        void insertSeg(SEG seg, short switch_level, linestringAttri attri);
        void insertFeature(LINESTRING line_g, short switch_level, LINESTRING_ITEM line);
        linestringNode *getMinWNode(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level);
        void countField_single(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, map<string, double> &count_map);
        void countField_classified(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, map<string, double> &count_map);
        void countField_ruled(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, vector<string> rules, map<string, double> &count_map);
        void getMaxFcount(int z, int &max_fcount);
};

class polygonNode{
    public:
        // node type
        nodeType node_type;
        char topo_type;
        // node BBox
        double node_x; 
        double node_y;
        double node_w;
        // node attributes
        unsigned short node_level;
        double node_farea; 
        // float ANN;
        // node ptrs
        POLYGON_ATTRI_NODE attri_map;
        // node ptrs
        polygonNode *parent;
        polygonNode *LUNode;
        polygonNode *RUNode;
        polygonNode *LBNode;
        polygonNode *RBNode;
        polygonRtree Rtree;
        // node functions
        polygonNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, polygonNode *_parent);
        ~polygonNode();
        polygonNode *insertTMBB(float obj_xMin, float obj_yMin, float obj_xMax, float obj_yMax, double farea, short switch_level, polygonAttri attri);
        void insertEdge(float edge_xMin, float edge_yMin, float edge_xMax, float edge_yMax, string edge_type, short switch_level, polygonAttri attri);
        void insertFeature(float box_xMin, float box_yMin, float box_xMax, float box_yMax, short switch_level, float ring_x[], float ring_y[], int len, polygonAttri attri, POLYGON poly, double farea, POLYGON_ITEM ring);
        polygonNode *getMinWNode(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level);
        void countField_single(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, map<string, double> &count_map);
        void countField_classified(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, map<string, double> &count_map);
        void countField_ruled(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, vector<string> rules, map<string, double> &count_map);
};


#endif