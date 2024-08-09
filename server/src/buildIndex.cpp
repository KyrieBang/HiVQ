#include "hicore/core.hpp"
#include "hicore/core_creator.hpp"
#include "hicore/core_iterator.hpp"
#include "hicore/shp.hpp"
#include "buildIndex.hpp"
#include "treeNode.hpp"
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <mpi.h>
#include <omp.h> 
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <regex>

using namespace std;
using namespace HiGIS::IO;
using namespace HiGIS::Core;

#define pi 3.14159265358979323
#define L 20037508.3427892

typedef pair<map<string, vector<string>>, map<string, vector<double>>> FIELDS_MAP;

// coordinate transformation (4326->3857)
void coordTran(double &x ,double &y){
    x = x * 20037508.34 / 180;
    y = (log(tan(((90 + y) * pi) / 360)) / (pi / 180)) * 20037508.34 / 180;
}


// calculate the spatial scale of the dataset
short calSpatialScale(GeoData data, string data_srs){
    short level = 0;
    double data_xmin = data.extent.xmin;
    double data_xmax = data.extent.xmax;
    double data_ymin = data.extent.ymin;
    double data_ymax = data.extent.ymax;
    if (data_srs == "4326"){
        coordTran(data_xmin, data_ymin);
        coordTran(data_xmax, data_ymax);
    }
    float dlt_x = data_xmax - data_xmin;
    float dlt_y = data_ymax - data_ymin;
    float dlt_max = max(dlt_x, dlt_y);
    float span = 2 * L;
    while (dlt_max < span){
        span = span / 2;
        level++;
    } 
    return level;
}


// get the info of field
void getFieldInfo(DBFHandle dbf, vector<int> &field_types, vector<int> &field_indexs, vector<string> &fields){
    char *field = new char[12];
    int *pnWidth, *pnDecimals;
    for (int i = 0; i < DBFGetFieldCount(dbf); i++){
        field_types.push_back(DBFGetFieldInfo(dbf, i, field, pnWidth, pnDecimals));
        fields.push_back((string) field);
        field_indexs.push_back(DBFGetFieldIndex(dbf, field));
    }
}


// get all fields of the vector dataset
template<typename T>
void getFields(int id, DBFHandle dbf, vector<int> field_types, vector<int> field_indexs, vector<string> fields, FieldSet p_fields, map<string, pair<int, int>> field_range, T &attri){
    char *field = new char[12];
    int *pnWidth, *pnDecimals;

    attri.id = id;

    for (auto field_c: p_fields.fields_classify){
        int ix = find(fields.begin(), fields.end(), field_c) - fields.begin();
        if (ix < fields.size()){
            string field_value = DBFReadStringAttribute(dbf, id, field_indexs[ix]);
            attri.field_string[field_c] = field_value;
        }
    }

    for (auto field_r: p_fields.fields_continuous){
        int ix = find(fields.begin(), fields.end(), field_r) - fields.begin();
        if (ix < fields.size()){
            int field_value = DBFReadIntegerAttribute(dbf, id, field_indexs[ix]);
            
            pair<int, int> min_max = field_range[field_r];
            int min_value = min_max.first;
            int max_value = min_max.second;
            int seg_num = 10;
            int seg_len = (max_value - min_value) / seg_num + 1;
            int class_value = (field_value - min_value) / seg_len;
            attri.field_string[field_r] = to_string(class_value * seg_len) + "~" + to_string((class_value + 1) * seg_len);
        }
    }



    // for (int i = 0; i < DBFGetFieldCount(dbf); i++){
    //     string field = fields[i];
    //     int ix = find(p_fields.begin(), p_fields.end(), field) - p_fields.begin();

    //     if (ix < p_fields.size()){
    //         switch(field_types[i]){
    //             case 0:
    //                 field_value = DBFReadStringAttribute(dbf, id, field_indexs[i]);
    //                 break;
    //             case 1:
    //                 field_value = to_string(DBFReadIntegerAttribute(dbf, id, field_indexs[i]));
    //                 break;
    //             case 2:
    //                 field_value = to_string(DBFReadDoubleAttribute(dbf, id, field_indexs[i]));
    //                 break;
    //         }

    //         attri.field_string[field] = field_value;
    //     }
    // }
}

map<string, pair<int, int>> getFieldMinMax(DBFHandle data_dbf, vector<int> field_indexs, vector<string> fields, FieldSet part_fields){
    map<string, pair<int, int>> field_range;
    vector<string> part_fields_c = part_fields.fields_continuous;
    for (auto field_real: part_fields_c){
        int ix = find(fields.begin(), fields.end(), field_real) - fields.begin();
        if (ix < fields.size()){
            vector<int> real_list;
            for (int i = 0; i < DBFGetRecordCount(data_dbf); i++){
                real_list.push_back(DBFReadIntegerAttribute(data_dbf, i, field_indexs[ix]));
            }
            int min_value = *min_element(real_list.begin(), real_list.end());
            int max_value = *max_element(real_list.begin(), real_list.end());
            field_range[field_real] = make_pair(min_value, max_value);
        }
    }
    return field_range;
}


// build index of point dataset
pair<pointNode*, short> pointIndex_shp(string data_path, string data_srs, FieldSet part_fields){
    struct timeval	t1, t2;
	gettimeofday(&t1, NULL);
    // init
    int data_count = 0;
    vector<int> field_indexs, field_types;
    vector<string> fields;
    // read shapfile
    Shapefile shp;
    GeoData data = shp.Read(data_path);
    DBFHandle data_dbf = DBFOpen(shp.AuxPath(data_path, ".dbf").c_str(), "r");
    getFieldInfo(data_dbf, field_types, field_indexs, fields);
    
    // get the min and max value of continuous attribution
    map<string, pair<int, int>> field_range = getFieldMinMax(data_dbf, field_indexs, fields, part_fields);
    
    // create root node ptr 
    pointNode *pointRNode = new pointNode(ROOT, -L, -L, 2 * L, 0, nullptr);
    // calculate the spatial scale of the dataset
    short data_level = calSpatialScale(data, data_srs);
    short switch_level = data_level + 9;
    // traverse point feature
    PointIterator ptIter(data);
    do{
        double x = ptIter.X();
        double y = ptIter.Y();
        if (data_srs == "4326"){
            coordTran(x, y);
        }

        pointAttri attri;
        getFields(data_count, data_dbf, field_types, field_indexs, fields, part_fields, field_range, attri);
        
        // insert point from root node
        pointRNode->insertObj(x, y, switch_level, attri);
        data_count++;
    }while (ptIter.NextFeature() >= 0);
    data_count = data.feature_count;
    gettimeofday(&t2, NULL);
    float time_use = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;

    DBFClose(data_dbf);

    // print information
    cout << "the number of point feature is: " << data_count << std::endl;
    cout << "the time of building index (" << switch_level << " + 8) is: " << time_use <<  " s" << endl;
    return make_pair(pointRNode, switch_level);
}


// build index of linestring dataset
pair<linestringNode*, short> linestringIndex_shp(string data_path, string data_srs, FieldSet part_fields){
    struct timeval	t1, t2;
	gettimeofday(&t1, NULL);
    // init
    int feature_count = 0;
    int segment_count = 0;
    float total_len = 0;
    vector<int> field_indexs, field_types;
    vector<string> fields;
    // read shapefile
    Shapefile shp;
    GeoData data = shp.Read(data_path);
    DBFHandle data_dbf = DBFOpen(shp.AuxPath(data_path, ".dbf").c_str(), "r");
    getFieldInfo(data_dbf, field_types, field_indexs, fields);
    // get the min and max value of continuous attribution
    map<string, pair<int, int>> field_range = getFieldMinMax(data_dbf, field_indexs, fields, part_fields);

    // create root node ptr
    linestringNode *lineRNode = new linestringNode(ROOT, -L, -L, 2 * L, 0, nullptr);
    // calculate the spatial scale of the dataset
    short data_level = calSpatialScale(data, data_srs);
    short switch_level = data_level + 9;
    // traverse linestring feature
    LineStringIterator lineIter(data);
    do{
        int nodes_num = lineIter.NodeCount();
        double line_x[nodes_num];
        double line_y[nodes_num];
        LINE line_coord;
        LINESTRING line_g;
        double *nodes = lineIter.Nodes();
        #pragma omp parallel for num_threads(4)
        for (int i = 0; i < nodes_num; i++){
            double xi = lineIter.X(nodes, i);
            double yi = lineIter.Y(nodes, i);
            if (data_srs == "4326"){
                coordTran(xi, yi);
            }
            line_x[i] = xi;
            line_y[i] = yi;
            line_coord.push_back(xi);
            line_coord.push_back(yi);

            bg::append(line_g, POINT(xi, yi));
        }
        double x_min = *min_element(line_x, line_x+nodes_num);
        double x_max = *max_element(line_x, line_x+nodes_num);
        double y_min = *min_element(line_y, line_y+nodes_num);
        double y_max = *max_element(line_y, line_y+nodes_num);
        // attri info
        linestringAttri attri;
        getFields(feature_count, data_dbf, field_types, field_indexs, fields, part_fields, field_range, attri);
        // insert TMBB from root node
        double line_length = bg::length(line_g)/1000.0;
        linestringNode *MBNode = lineRNode->insertTMBB(x_min, y_min, x_max, y_max, (nodes_num-1), line_length, switch_level, attri);
        // insert each segement from the minimum bounding node
        for (int j = 0; j < nodes_num - 1; j++){
            double x0 = line_x[j];
            double y0 = line_y[j];
            double x1 = line_x[j + 1];
            double y1 = line_y[j + 1];
            MBNode->insertSeg(SEG(POINT(x0, y0), POINT(x1, y1)), switch_level, attri);
        }
        // inesrt feature
        LINESTRING_ITEM line = make_tuple(Box(POINT(x_min, y_min), POINT(x_max, y_max)), line_coord, attri);
        MBNode->insertFeature(line_g, switch_level, line);

        segment_count += (nodes_num - 1);
        feature_count++;
        total_len += line_length;
    }while(lineIter.NextFeature() >= 0);
    feature_count = data.feature_count;
    gettimeofday(&t2, NULL);
    float time_use = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;

    DBFClose(data_dbf);

    // print information
    cout << "the number of linestring feature is: " << feature_count << endl;
    cout << "the number of the total segement is: " << segment_count << endl;
    cout << "the total length of linestring feature is: " << total_len << " KM" << endl;
    cout << "the time of building index (" << switch_level << " + 9) is: " << time_use <<  " s" << endl;

    return make_pair(lineRNode, switch_level);
}


// build index of polygon dataset
pair<polygonNode*, short> polygonIndex_shp(string data_path, string data_srs, FieldSet part_fields){
    struct timeval	t1, t2;
	gettimeofday(&t1, NULL);
    // init
    int ring_count = 0;
    int feature_count = 0;
    int segment_count = 0;
    double total_area = 0;
    vector<int> field_indexs, field_types;
    vector<string> fields;
    // read shapefile
    Shapefile shp;
    GeoData data = shp.Read(data_path);
    auto data_dbf = DBFOpen(shp.AuxPath(data_path, ".dbf").c_str(), "r");
    getFieldInfo(data_dbf, field_types, field_indexs, fields);
    // get the min and max value of continuous attribution
    map<string, pair<int, int>> field_range = getFieldMinMax(data_dbf, field_indexs, fields, part_fields);

    // create root node
    polygonNode *polygonRNode = new polygonNode(ROOT, -L, -L, 2 * L, 0, nullptr);
    // calculate the spatial scale of the dataset
    short data_level = calSpatialScale(data, data_srs);
    // short switch_level = data_level + 9;
    short switch_level = 12;
    // traverse polygon feature
    PolygonIterator polygonIter(data);
    do{
        // get the polygon feature
        for (int n = 0; n < polygonIter.PartCount(); n++){
            auto poly = polygonIter.Part(n);
            // get the ring of each polygon
            for (int k = 0; k < polygonIter.RingCount(poly); k++){
                auto ring = polygonIter.Ring(poly, k);
                auto nodes = polygonIter.Nodes(ring);
                int nodes_num = polygonIter.NodeCount(ring);
                float ring_x[nodes_num];
                float ring_y[nodes_num];
                RING ring_cord;
                POLYGON polygon_g;
                // for each node of each ring
                for (int i = 0; i < nodes_num; i++){
                    double xi = polygonIter.X(nodes, i);
                    double yi = polygonIter.Y(nodes, i);
                    if (data_srs == "4326"){
                        coordTran(xi, yi);
                    }
                    ring_x[i] = xi;
                    ring_y[i] = yi;
                    ring_cord.push_back(xi);
                    ring_cord.push_back(yi);

                    polygon_g.outer().push_back(POINT(xi, yi));
                }
                double x_min = *min_element(ring_x, ring_x+nodes_num);
                double x_max = *max_element(ring_x, ring_x+nodes_num);
                double y_min = *min_element(ring_y, ring_y+nodes_num);
                double y_max = *max_element(ring_y, ring_y+nodes_num);
                // init field struct
                polygonAttri attri;
                getFields(feature_count, data_dbf, field_types, field_indexs, fields, part_fields, field_range, attri);
                // insert TMBB from root node
                double farea = abs(bg::area(polygon_g) / 1000000);
                total_area += farea;
                polygonNode *MBNode = polygonRNode->insertTMBB(x_min, y_min, x_max, y_max, farea, switch_level, attri);
                // insert each edge from the minimum bounding node
                for (int j = 0; j < nodes_num - 1; j++){
                    double x0 = ring_x[j];
                    double y0 = ring_y[j];
                    double x1 = ring_x[j + 1];
                    double y1 = ring_y[j + 1];
                    float edge_xMin = min(x0, x1);
                    float edge_xMax = max(x0, x1);
                    float edge_yMin = min(y0, y1);
                    float edge_yMax = max(y0, y1);
                    // POLYGON poly = make_pair(SEG(POINT(x0, y0), POINT(x1, y1)), ring_cord); 
                    if (((x1 >= x0) && (y1 >= y0)) or ((x1 <= x0) && (y1 <= y0))){
                        MBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, "0", switch_level, attri);
                    }
                    else if (((x1 < x0) && (y1 > y0)) or ((x1 > x0) && (y1 < y0))){
                        MBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, "1", switch_level, attri);
                    }
                }
                
                // insert the polygon feature from the minimum bounding node
                POLYGON_ITEM poly_ring = make_tuple(Box(POINT(x_min, y_min), POINT(x_max, y_max)), ring_cord, attri);
                MBNode->insertFeature(x_min, y_min, x_max, y_max, switch_level, ring_x, ring_y, nodes_num, attri, polygon_g, farea, poly_ring);
                segment_count += (nodes_num-1);
                // ring_count++;
            }
        }
        feature_count++;
    }while(polygonIter.NextFeature() >= 0);
    feature_count = data.feature_count;
    gettimeofday(&t2, NULL);
    float time_use = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;

    DBFClose(data_dbf);

    // print information
    cout << "the number of polyogn feature is: " << feature_count << endl;
    cout << "the number of the total edge is: " << segment_count << endl;
    cout << "the number of the total area is: " << total_area << " KM2" << endl;
    cout << "the size of building index (" << switch_level << " + 9) is: " << time_use <<  " s" << endl;
    return make_pair(polygonRNode, switch_level);
}

