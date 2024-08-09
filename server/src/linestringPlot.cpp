#include <iostream>
#include <fstream>
#include "png.h"
#include "linestringPlot.hpp"
#include "treeNode.hpp"
#include "pointPlot.hpp"
#include "intersectJudge.hpp"
#include <omp.h>
#include <mpi.h>
#include <mapnik/map.hpp>
#include <mapnik/agg_renderer.hpp>
#include <mapnik/image.hpp>
#include <mapnik/image_util.hpp>
#include <mapnik/layer.hpp>
#include <mapnik/rule.hpp>
#include <mapnik/feature_type_style.hpp>
#include <mapnik/symbolizer.hpp>
#include <mapnik/text/placements/dummy.hpp>
#include <mapnik/text/text_properties.hpp>
#include <mapnik/text/formatting/text.hpp>
#include <mapnik/datasource_cache.hpp>
#include <mapnik/datasource.hpp>
#include <mapnik/font_engine_freetype.hpp>
#include <mapnik/expression.hpp>
#include <mapnik/color_factory.hpp>
#include <mapnik/unicode.hpp>
#include <mapnik/save_map.hpp>
#include <mapnik/cairo_io.hpp>
#include <mapnik/geometry.hpp>
#include <mapnik/geometry_envelope.hpp>
#include <mapnik/query.hpp>
#include <mapnik/geometry_centroid.hpp>
using namespace mapnik;
using namespace mapnik::geometry;
#define L 20037508.3427892
#define TILE_SIZE 256


const bool mapnik_statue = datasource_cache::instance().register_datasources("/usr/lib/mapnik/3.0/input");
const string srs_merc = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0.0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs +over";


bool cmp_value(const pair<string, double> left, const pair<string, double> right)
{
	return left.second < right.second;
}

int attriClassify(map<string, double> node_attri, vector<string> attri_list){
    map<string, double> tmp;
    for (int i = 0; i < attri_list.size(); i++){
        string attri_i = attri_list[i];
        if (node_attri.find(attri_i) != node_attri.end()){
            tmp[attri_i] = node_attri[attri_i];
        }
    }

    if (!tmp.size()){
        return attri_list.size();
    }
    else {
        auto max_pair = max_element(tmp.begin(), tmp.end(), cmp_value);
        return (find(attri_list.begin(), attri_list.end(), max_pair->first)-attri_list.begin());
    }
}

int attriRule(map<string, double> node_attri, vector<string> filter_list){
    string attri;
    map<string, double> tmp;
    map<string, int> tmp2;
    for (int i = 0; i < filter_list.size(); i++){
        string filter_i = filter_list[i];
        int ix;
        for (auto node_i: node_attri){
            string node_i_key = node_i.first;
            // int node_i_value = node_i.second.first;
            if ((ix = filter_i.find("!=")) != string::npos){
                if (node_i_key != filter_i.substr(ix+2)){
                    tmp[node_i_key] = node_i.second;
                    tmp2[node_i_key] = i;
                }
            }
            else if ((ix = filter_i.find("%3E=")) != string::npos){
                if (stod(node_i_key) >= stod(filter_i.substr(ix+4))){
                    tmp[node_i_key] = node_i.second;
                    tmp2[node_i_key] = i;
                }
            }
            else if ((ix = filter_i.find("%3C=")) != string::npos){
                if (stod(node_i_key) <= stod(filter_i.substr(ix+4))){
                    tmp[node_i_key] = node_i.second;
                    tmp2[node_i_key] = i;
                }
            }
            else if ((ix = filter_i.find('=')) != string::npos){
                if (node_i_key == filter_i.substr(ix+1)){
                    tmp[node_i_key] = node_i.second;
                    tmp2[node_i_key] = i;
                }
            }
            else if ((ix = filter_i.find("%3E")) != string::npos){
                if (stod(node_i_key) > stod(filter_i.substr(ix+3))){
                    tmp[node_i_key] = node_i.second;
                    tmp2[node_i_key] = i;
                }
            }
            else if ((ix = filter_i.find("%3C")) != string::npos){
                if (stod(node_i_key) < stod(filter_i.substr(ix+3))){
                    tmp[node_i_key] = node_i.second;
                    tmp2[node_i_key] = i;
                }
            }
        }
    }

    if (!tmp.size()){
        return filter_list.size();
    }
    else{
        auto max_pair = max_element(tmp.begin(), tmp.end(), cmp_value);
        attri = max_pair->first;
    }
    
    return tmp2[attri];
}


// determine whether the node correspongding to the pixel exists
bool nodeIfExist(int z, double tile_box[], double pix_x, double pix_y, linestringNode *seek){
    double min_x = tile_box[0];
    double min_y = tile_box[1];
    double max_x = tile_box[2];
    double max_y = tile_box[3];
    double mid_x, mid_y;
    int divided_num = 0;

    while (divided_num < 8){
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // pixel in LU sub-region
        if (pix_x < mid_x && pix_y >= mid_y){
            if (seek->LUNode == nullptr){
                return false;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // pixel in RU sub-region
        else if (pix_x >= mid_x && pix_y >= mid_y){
            if (seek->RUNode == nullptr){
                return false;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // pixel in LB sub-region
        else if (pix_x < mid_x && pix_y < mid_y){
            if (seek->LBNode == nullptr){
                return false;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // pixel in RB sub-region
        else if (pix_x >= mid_x && pix_y < mid_y){
            if (seek->RBNode == nullptr){
                return false;
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }else {
            return false;
        }
        divided_num ++;
    }
    return true;
}


linestringNode * nodeIfExist_attri(int z, double tile_box[], double pix_x, double pix_y, linestringNode *seek){
    double min_x = tile_box[0];
    double min_y = tile_box[1];
    double max_x = tile_box[2];
    double max_y = tile_box[3];
    double mid_x, mid_y;
    int divided_num = 0;

    while (divided_num < 8){
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // pixel in LU sub-region
        if (pix_x < mid_x && pix_y >= mid_y){
            if (seek->LUNode == nullptr){
                return nullptr;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // pixel in RU sub-region
        else if (pix_x >= mid_x && pix_y >= mid_y){
            if (seek->RUNode == nullptr){
                return nullptr;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // pixel in LB sub-region
        else if (pix_x < mid_x && pix_y < mid_y){
            if (seek->LBNode == nullptr){
                return nullptr;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // pixel in RB sub-region
        else if (pix_x >= mid_x && pix_y < mid_y){
            if (seek->RBNode == nullptr){
                return nullptr;
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }else {
            return nullptr;
        }
        divided_num ++;
    }

    return seek;
}



// determine whether the node correspongding to the tile<z,x,y> exists
linestringNode * locTileNode(int z, int x, int y, linestringNode *seek){
    double tile_Rz = L / pow(2, z-1);
    double tile_x = (x + 0.5) * tile_Rz - L;
    double tile_y = L - (y + 0.5) * tile_Rz;
    double min_x = -L;
    double max_x = L;
    double min_y = -L;
    double max_y = L;
    double mid_x, mid_y;
    int divided_num = 0;

    while (divided_num < z){
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // tile_center_pt in LU sub-region
        if (tile_x < mid_x && tile_y >= mid_y){
            if (seek->LUNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in RU sub-region
        else if (tile_x >= mid_x && tile_y >= mid_y){
            if (seek->RUNode == nullptr){
                seek = nullptr;
                break;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in LB sub-region
        else if (tile_x < mid_x && tile_y < mid_y){
            if (seek->LBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // tile_center_pt in RB sub-region
        else if (tile_x >= mid_x && tile_y < mid_y){
            if (seek->RBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }
        divided_num ++;
    }
    return seek;
}


// find the node which links to rtree
linestringNode * locRtreeNode(int did_num, double tile_x, double tile_y, linestringNode *seek, short switch_level){
    double min_x = -L;
    double max_x = L;
    double min_y = -L;
    double max_y = L;
    double mid_x, mid_y;
    int divided_num = 0;
    linestringNode *middle_node = nullptr;
    if (did_num > switch_level + 9)
        did_num = switch_level + 9;

    while (divided_num < did_num){
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // tile_center_pt in LU sub-region
        if (tile_x < mid_x && tile_y >= mid_y){
            if (seek->LUNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in RU sub-region
        else if (tile_x >= mid_x && tile_y >= mid_y){
            if (seek->RUNode == nullptr){
                seek = nullptr;
                break;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in LB sub-region
        else if (tile_x < mid_x && tile_y < mid_y){
            if (seek->LBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // tile_center_pt in RB sub-region
        else if (tile_x >= mid_x && tile_y < mid_y){
            if (seek->RBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }
        
        // get the node which links to r-tree
        if (divided_num == switch_level){
            middle_node = seek;
        }

        divided_num ++;
    }

    if (seek){
        return middle_node;
    }
    return seek;
}


// determine whether the node correspongding to the nearest pixel exists
linestringNode * locPixNode(int z, double tile_box[], double pix_x, double pix_y, linestringNode *seek, linestringNode *RNode){
    double min_x = tile_box[0];
    double min_y = tile_box[1];
    double max_x = tile_box[2];
    double max_y = tile_box[3];
    double mid_x, mid_y;

    if (abs(pix_x) <= L && abs(pix_y) <= L){
        if (pix_x < min_x || pix_x > max_x || pix_y < min_y || pix_y > max_y){
            seek = RNode;
            min_x = -L;
            max_x = L;
            min_y = -L;
            max_y = L;
        }
    }
    else{
        seek = nullptr;
        return seek;
    }

    int level = 8 + z - seek->node_level;
    int divided_num = 0;

    while (divided_num < level){
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // pixel in LU sub-region
        if (pix_x < mid_x && pix_y >= mid_y){
            if (seek->LUNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // pixel in RU sub-region
        else if (pix_x >= mid_x && pix_y >= mid_y){
            if (seek->RUNode == nullptr){
                seek = nullptr;
                break;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // pixel in LB sub-region
        else if (pix_x < mid_x && pix_y < mid_y){
            if (seek->LBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // pixel in RB sub-region
        else if (pix_x >= mid_x && pix_y < mid_y){
            if (seek->RBNode == nullptr){
                seek = nullptr;
                break;
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }
        divided_num++;
    }
    return seek;
}


// get the child node of the node which is nearest to the pixel 
double *getNearestCNode(linestringNode *seek, double pix_x, double pix_y, double near_x, double near_y, double Rz){
    double dist = 0;
    double node_x, node_y;
    if (seek->LUNode != nullptr){
        double node_x1 = near_x - Rz;
        double node_y1 = near_y + Rz;
        double dist1 = sqrt(pow(node_x1 - pix_x, 2) + pow(node_y1 - pix_y, 2));
        if (dist == 0){
            dist = dist1;
            node_x = node_x1;
            node_y = node_y1;
        }
    }
    if (seek->RUNode != nullptr){
        double node_x2 = near_x + Rz;
        double node_y2 = near_y + Rz;
        double dist2 = sqrt(pow(node_x2 - pix_x, 2) + pow(node_y2 - pix_y, 2));
        if (dist == 0 || dist2 < dist){
            dist = dist2;
            node_x = node_x2;
            node_y = node_y2;
        }
    }
    if (seek->LBNode != nullptr){
        double node_x3 = near_x - Rz;
        double node_y3 = near_y - Rz;
        double dist3 = sqrt(pow(node_x3 - pix_x, 2) + pow(node_y3 - pix_y, 2));
        if (dist == 0 || dist3 < dist){
            dist = dist3;
            node_x = node_x3;
            node_y = node_y3;
        }
    }
    if (seek->RBNode != nullptr){
        double node_x4 = near_x + Rz;
        double node_y4 = near_y - Rz;
        double dist4 = sqrt(pow(node_x4 - pix_x, 2) + pow(node_y4 - pix_y, 2));
        if (dist == 0 || dist4 < dist){
            dist = dist4;
            node_x = node_x4;
            node_y = node_y4;
        }
    }
    double *nearInfo = new double[3];
    nearInfo[0] = dist;
    nearInfo[1] = node_x;
    nearInfo[2] = node_y;
    return nearInfo;
}


// display-driven tile plotting
void dis_tileRender_single(int z, int x, int y, int R, int G, int B, float AD, float width, linestringNode *tile_node, linestringNode *linestringRNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    double R0 = width * Rz; 
    double R1 = R0 - sqrt(2) * Rz / 4.0;
    double R2 = R0 + sqrt(2) * Rz / 4.0;
    double Rz4 = 0.25 * Rz;
    int level = (width == 0) ? 0 : int(width + sqrt(2) / 4.0 - 0.5) + 1;
    // calcalate the BBox of tile <z,x,y>
    double tile_Rz = 256 * Rz;
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            // determine whether the node correspongding to the pixel exists
            if (nodeIfExist(z, tile_box, pix_x, pix_y, tile_node)){
                png_ptr[i][4 * j] = R;
                png_ptr[i][4 * j + 1] = G;
                png_ptr[i][4 * j + 2] = B;
                png_ptr[i][4 * j + 3] = AD * 255;
                continue;
            }

            short final_value = 0;
            short pix_value = 0;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (m == -l || m == l || n == -l || n == l){
                            linestringNode *node_floor = locPixNode(z, tile_box, pix_x + m * Rz, pix_y + n * Rz, tile_node, linestringRNode);
                            if (node_floor != nullptr){
                                double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, pix_x + m * Rz, pix_y + n * Rz, Rz / 4);
                                double dist = nearInfo[0];
                                if (dist < R1){
                                    final_value = 4;
                                    break;
                                }
                                else if (dist < R2){
                                    double node_x = nearInfo[1];
                                    double node_y = nearInfo[2];
                                    if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0){
                                        pix_value += 1;
                                    }
                                    if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0){
                                        pix_value += 1;
                                    }
                                    if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0){
                                        pix_value += 1;
                                    }  
                                    if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0){
                                        pix_value += 1; 
                                    }

                                    if (pix_value == 4){
                                        final_value = 4;
                                        break;
                                    }
                                    else if (pix_value > final_value){
                                        final_value = pix_value;
                                    }
                                    pix_value = 0; 
                                }
                            }
                        }
                        
                    }
                    if (final_value == 4)
                        break;
                }
                if (final_value == 4)
                    break;
            }
            // generate pixel value
            if (final_value > 0){
                png_ptr[i][4 * j] = R;
                png_ptr[i][4 * j + 1] = G;
                png_ptr[i][4 * j + 2] = B;
                png_ptr[i][4 * j + 3] = AD * final_value * 64 - 1;
            }
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
        }
    }
}

void dis_tileRender_range_single(int z, int x, int y, double box[], int R, int G, int B, float AD, float width, linestringNode *tile_node, linestringNode *linestringRNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    double R0 = width * Rz; 
    double R1 = R0 - sqrt(2) * Rz / 4.0;
    double R2 = R0 + sqrt(2) * Rz / 4.0;
    double Rz4 = 0.25 * Rz;
    int level = (width == 0) ? 0 : int(width + sqrt(2) / 4.0 - 0.5) + 1;
    // calcalate the BBox of tile <z,x,y>
    double tile_Rz = 256 * Rz;
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    double box_xMin = box[0];
    double box_yMin = box[1];
    double box_xMax = box[2];
    double box_yMax = box[3];
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_xMin = (256 * x + j) * Rz - L;
            double pix_xMax = (256 * x + j + 1) * Rz - L;
            double pix_yMin = L - (256 * y + i + 1) * Rz;
            double pix_yMax = L - (256 * y + i) * Rz;
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            if (IfOverlapMBB(pix_xMin, pix_yMin, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                // determine whether the node correspongding to the pixel exists
                if (nodeIfExist(z, tile_box, pix_x, pix_y, tile_node)){
                    png_ptr[i][4 * j] = R;
                    png_ptr[i][4 * j + 1] = G;
                    png_ptr[i][4 * j + 2] = B;
                    png_ptr[i][4 * j + 3] = AD * 255;
                    continue;
                }
            }

            short final_value = 0;
            short pix_value = 0;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (m == -l || m == l || n == -l || n == l){
                            if (IfOverlapMBB(pix_xMin + m * Rz, pix_yMin + n * Rz, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                                linestringNode *node_floor = locPixNode(z, tile_box, pix_x + m * Rz, pix_y + n * Rz, tile_node, linestringRNode);
                                if (node_floor != nullptr){
                                    double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, pix_x + m * Rz, pix_y + n * Rz, Rz / 4);
                                    double dist = nearInfo[0];
                                    if (dist < R1){
                                        final_value = 4;
                                        break;
                                    }
                                    else if (dist < R2){
                                        double node_x = nearInfo[1];
                                        double node_y = nearInfo[2];
                                        if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0){
                                            pix_value += 1;
                                        }
                                        if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0){
                                            pix_value += 1;
                                        }
                                        if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0){
                                            pix_value += 1;
                                        }  
                                        if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0){
                                            pix_value += 1; 
                                        }

                                        if (pix_value == 4){
                                            final_value = 4;
                                            break;
                                        }
                                        else if (pix_value > final_value){
                                            final_value = pix_value;
                                        }
                                        pix_value = 0; 
                                    }
                                }
                            }
                        }
                        
                    }
                    if (final_value == 4)
                        break;
                }
                if (final_value == 4)
                    break;
            }
            // generate pixel value
            if (final_value > 0){
                png_ptr[i][4 * j] = R;
                png_ptr[i][4 * j + 1] = G;
                png_ptr[i][4 * j + 2] = B;
                png_ptr[i][4 * j + 3] = AD * final_value * 64 - 1;
            }
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
        }
    }
}

void dis_tileRender_classify(int z, int x, int y, vector<string> attri, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string field_name, linestringNode *tile_node, linestringNode *linestringRNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    vector<double> R0, R1, R2;
    for (int k = 0; k < width.size(); k++){
        R0.push_back(width[k] * Rz);
        R1.push_back(width[k] * Rz - sqrt(2) * Rz / 4.0);
        R2.push_back(width[k] * Rz + sqrt(2) * Rz / 4.0);
    }
    double Rz4 = 0.25 * Rz;
    float max_width = *max_element(width.begin(), width.end());
    int level = (max_width == 0) ? 0 : int(max_width + sqrt(2) / 4.0 - 0.5) + 1;
    // calcalate the BBox of tile <z,x,y>
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            // determine whether the node correspongding to the pixel exists
            linestringNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
            if (pixel_node != nullptr){
                int ix = attriClassify(pixel_node->attri_map[field_name], attri);
                if (ix != attri.size()){
                    png_ptr[i][4 * j] = R[ix];     
                    png_ptr[i][4 * j + 1] = G[ix]; 
                    png_ptr[i][4 * j + 2] = B[ix]; 
                    png_ptr[i][4 * j + 3] = AD[ix] * 255;
                    continue;
                }
            }
            
            short final_value = 0;
            int final_ix;
            short pix_value = 0;
            int attri_ix;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (m == -l || m == l || n == -l || n == l){
                            linestringNode *node_floor = locPixNode(z, tile_box, pix_x + m * Rz, pix_y + n * Rz, tile_node, linestringRNode);
                            if (node_floor != nullptr){
                                attri_ix = attriClassify(node_floor->attri_map[field_name], attri);
                                if (attri_ix == attri.size()){
                                    continue;
                                }

                                double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, pix_x + m * Rz, pix_y + n * Rz, Rz / 4);
                                double dist = nearInfo[0];
                                if (dist < R1[attri_ix]){
                                    final_value = 4;
                                    final_ix = attri_ix;
                                    break;
                                }
                                else if (dist < R2[attri_ix]){
                                    double node_x = nearInfo[1];
                                    double node_y = nearInfo[2];
                                    if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
                                        pix_value += 1;
                                    }
                                    if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
                                        pix_value += 1;
                                    }
                                    if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
                                        pix_value += 1;
                                    }  
                                    if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
                                        pix_value += 1; 
                                    }

                                    if (pix_value == 4){
                                        final_value = 4;
                                        final_ix = attri_ix;
                                        break;
                                    }
                                    else if (pix_value > final_value){
                                        final_value = pix_value;
                                        final_ix = attri_ix;
                                    }
                                    pix_value = 0; 
                                }
                            }
                        }
                        
                    }
                    if (final_value == 4)
                        break;
                }
                if (final_value == 4)
                    break;
            }
            // generate pixel value
            if (final_value > 0 and final_ix != attri.size()){
                png_ptr[i][4 * j] = R[final_ix];
                png_ptr[i][4 * j + 1] = G[final_ix];
                png_ptr[i][4 * j + 2] = B[final_ix];
                png_ptr[i][4 * j + 3] = AD[final_ix] * final_value * 64 - 1;
            }
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
        }
    }
}

void dis_tileRender_range_classify(int z, int x, int y, double box[], vector<string> attri, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string field_name, linestringNode *tile_node, linestringNode *linestringRNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    vector<double> R0, R1, R2;
    for (int k = 0; k < width.size(); k++){
        R0.push_back(width[k] * Rz);
        R1.push_back(width[k] * Rz - sqrt(2) * Rz / 4.0);
        R2.push_back(width[k] * Rz + sqrt(2) * Rz / 4.0);
    }
    double Rz4 = 0.25 * Rz;
    float max_width = *max_element(width.begin(), width.end());
    int level = (max_width == 0) ? 0 : int(max_width + sqrt(2) / 4.0 - 0.5) + 1;
    // calcalate the BBox of tile <z,x,y>
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    double box_xMin = box[0];
    double box_yMin = box[1];
    double box_xMax = box[2];
    double box_yMax = box[3];
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_xMin = (256 * x + j) * Rz - L;
            double pix_xMax = (256 * x + j + 1) * Rz - L;
            double pix_yMin = L - (256 * y + i + 1) * Rz;
            double pix_yMax = L - (256 * y + i) * Rz;
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            if (IfOverlapMBB(pix_xMin, pix_yMin, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                // determine whether the node correspongding to the pixel exists
                linestringNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
                if (pixel_node != nullptr){
                    int ix = attriClassify(pixel_node->attri_map[field_name], attri);
                    if (ix != attri.size()){
                        png_ptr[i][4 * j] = R[ix];     
                        png_ptr[i][4 * j + 1] = G[ix]; 
                        png_ptr[i][4 * j + 2] = B[ix]; 
                        png_ptr[i][4 * j + 3] = AD[ix] * 255;
                        continue;
                    }
                }
            }
            
            short final_value = 0;
            int final_ix;
            short pix_value = 0;
            int attri_ix;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (m == -l || m == l || n == -l || n == l){
                            if (IfOverlapMBB(pix_xMin + m * Rz, pix_yMin + n * Rz, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                                linestringNode *node_floor = locPixNode(z, tile_box, pix_x + m * Rz, pix_y + n * Rz, tile_node, linestringRNode);
                                if (node_floor != nullptr){
                                    attri_ix = attriClassify(node_floor->attri_map[field_name], attri);
                                    if (attri_ix == attri.size()){
                                        continue;
                                    }

                                    double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, pix_x + m * Rz, pix_y + n * Rz, Rz / 4);
                                    double dist = nearInfo[0];
                                    if (dist < R1[attri_ix]){
                                        final_value = 4;
                                        final_ix = attri_ix;
                                        break;
                                    }
                                    else if (dist < R2[attri_ix]){
                                        double node_x = nearInfo[1];
                                        double node_y = nearInfo[2];
                                        if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
                                            pix_value += 1;
                                        }
                                        if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
                                            pix_value += 1;
                                        }
                                        if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
                                            pix_value += 1;
                                        }  
                                        if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
                                            pix_value += 1; 
                                        }

                                        if (pix_value == 4){
                                            final_value = 4;
                                            final_ix = attri_ix;
                                            break;
                                        }
                                        else if (pix_value > final_value){
                                            final_value = pix_value;
                                            final_ix = attri_ix;
                                        }
                                        pix_value = 0; 
                                    }
                                }
                            }
                        }
                        
                    }
                    if (final_value == 4)
                        break;
                }
                if (final_value == 4)
                    break;
            }
            // generate pixel value
            if (final_value > 0 and final_ix != attri.size()){
                png_ptr[i][4 * j] = R[final_ix];
                png_ptr[i][4 * j + 1] = G[final_ix];
                png_ptr[i][4 * j + 2] = B[final_ix];
                png_ptr[i][4 * j + 3] = AD[final_ix] * final_value * 64 - 1;
            }
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
        }
    }
}

void dis_tileRender_rule(int z, int x, int y, vector<string> filters, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string field_name, linestringNode *tile_node, linestringNode *linestringRNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    vector<double> R0, R1, R2;
    for (int k = 0; k < width.size(); k++){
        R0.push_back(width[k] * Rz);
        R1.push_back(width[k] * Rz - sqrt(2) * Rz / 4.0);
        R2.push_back(width[k] * Rz + sqrt(2) * Rz / 4.0);
    }
    double Rz4 = 0.25 * Rz;
    float max_width = *max_element(width.begin(), width.end());
    int level = (max_width == 0) ? 0 : int(max_width + sqrt(2) / 4.0 - 0.5) + 1;
    // calcalate the BBox of tile <z,x,y>
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            // determine whether the node correspongding to the pixel exists
            linestringNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
            if (pixel_node != nullptr){
                int ix = attriRule(pixel_node->attri_map[field_name], filters);
                if (ix != filters.size()){
                    png_ptr[i][4 * j] = R[ix];     
                    png_ptr[i][4 * j + 1] = G[ix]; 
                    png_ptr[i][4 * j + 2] = B[ix]; 
                    png_ptr[i][4 * j + 3] = AD[ix] * 255;
                    continue;
                }
            }
            
            short final_value = 0;
            int final_ix;
            short pix_value = 0;
            int attri_ix;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (m == -l || m == l || n == -l || n == l){
                            linestringNode *node_floor = locPixNode(z, tile_box, pix_x + m * Rz, pix_y + n * Rz, tile_node, linestringRNode);
                            if (node_floor != nullptr){
                                attri_ix = attriRule(node_floor->attri_map[field_name], filters);
                                if (attri_ix == filters.size()){
                                    continue;
                                }

                                double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, pix_x + m * Rz, pix_y + n * Rz, Rz / 4);
                                double dist = nearInfo[0];
                                if (dist < R1[attri_ix]){
                                    final_value = 4;
                                    final_ix = attri_ix;
                                    break;
                                }
                                else if (dist < R2[attri_ix]){
                                    double node_x = nearInfo[1];
                                    double node_y = nearInfo[2];
                                    if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
                                        pix_value += 1;
                                    }
                                    if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
                                        pix_value += 1;
                                    }
                                    if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
                                        pix_value += 1;
                                    }  
                                    if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
                                        pix_value += 1; 
                                    }

                                    if (pix_value == 4){
                                        final_value = 4;
                                        final_ix = attri_ix;
                                        break;
                                    }
                                    else if (pix_value > final_value){
                                        final_value = pix_value;
                                        final_ix = attri_ix;
                                    }
                                    pix_value = 0; 
                                }
                            }
                        }
                        
                    }
                    if (final_value == 4)
                        break;
                }
                if (final_value == 4)
                    break;
            }
            // generate pixel value
            if (final_value > 0 and final_ix != filters.size()){
                png_ptr[i][4 * j] = R[final_ix];
                png_ptr[i][4 * j + 1] = G[final_ix];
                png_ptr[i][4 * j + 2] = B[final_ix];
                png_ptr[i][4 * j + 3] = AD[final_ix] * final_value * 64 - 1;
            }
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
        }
    }
}

void dis_tileRender_range_rule(int z, int x, int y, double box[], vector<string> filters, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string field_name, linestringNode *tile_node, linestringNode *linestringRNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    vector<double> R0, R1, R2;
    for (int k = 0; k < width.size(); k++){
        R0.push_back(width[k] * Rz);
        R1.push_back(width[k] * Rz - sqrt(2) * Rz / 4.0);
        R2.push_back(width[k] * Rz + sqrt(2) * Rz / 4.0);
    }
    double Rz4 = 0.25 * Rz;
    float max_width = *max_element(width.begin(), width.end());
    int level = (max_width == 0) ? 0 : int(max_width + sqrt(2) / 4.0 - 0.5) + 1;
    // calcalate the BBox of tile <z,x,y>
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    double box_xMin = box[0];
    double box_yMin = box[1];
    double box_xMax = box[2];
    double box_yMax = box[3];
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_xMin = (256 * x + j) * Rz - L;
            double pix_xMax = (256 * x + j + 1) * Rz - L;
            double pix_yMin = L - (256 * y + i + 1) * Rz;
            double pix_yMax = L - (256 * y + i) * Rz;
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            if (IfOverlapMBB(pix_xMin, pix_yMin, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                // determine whether the node correspongding to the pixel exists
                linestringNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
                if (pixel_node != nullptr){
                    int ix = attriRule(pixel_node->attri_map[field_name], filters);
                    if (ix != filters.size()){
                        png_ptr[i][4 * j] = R[ix];     
                        png_ptr[i][4 * j + 1] = G[ix]; 
                        png_ptr[i][4 * j + 2] = B[ix]; 
                        png_ptr[i][4 * j + 3] = AD[ix] * 255;
                        continue;
                    }
                }
            }
            
            short final_value = 0;
            int final_ix;
            short pix_value = 0;
            int attri_ix;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (m == -l || m == l || n == -l || n == l){
                            if (IfOverlapMBB(pix_xMin + m * Rz, pix_yMin + n * Rz, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                                linestringNode *node_floor = locPixNode(z, tile_box, pix_x + m * Rz, pix_y + n * Rz, tile_node, linestringRNode);
                                if (node_floor != nullptr){
                                    attri_ix = attriRule(node_floor->attri_map[field_name], filters);
                                    if (attri_ix == filters.size()){
                                        continue;
                                    }

                                    double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, pix_x + m * Rz, pix_y + n * Rz, Rz / 4);
                                    double dist = nearInfo[0];
                                    if (dist < R1[attri_ix]){
                                        final_value = 4;
                                        final_ix = attri_ix;
                                        break;
                                    }
                                    else if (dist < R2[attri_ix]){
                                        double node_x = nearInfo[1];
                                        double node_y = nearInfo[2];
                                        if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
                                            pix_value += 1;
                                        }
                                        if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
                                            pix_value += 1;
                                        }
                                        if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
                                            pix_value += 1;
                                        }  
                                        if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
                                            pix_value += 1; 
                                        }

                                        if (pix_value == 4){
                                            final_value = 4;
                                            final_ix = attri_ix;
                                            break;
                                        }
                                        else if (pix_value > final_value){
                                            final_value = pix_value;
                                            final_ix = attri_ix;
                                        }
                                        pix_value = 0; 
                                    }
                                }
                            }
                        }
                        
                    }
                    if (final_value == 4)
                        break;
                }
                if (final_value == 4)
                    break;
            }
            // generate pixel value
            if (final_value > 0 and final_ix != filters.size()){
                png_ptr[i][4 * j] = R[final_ix];
                png_ptr[i][4 * j + 1] = G[final_ix];
                png_ptr[i][4 * j + 2] = B[final_ix];
                png_ptr[i][4 * j + 3] = AD[final_ix] * final_value * 64 - 1;
            }
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
        }
    }
}

void dis_tileRender_heat(int z, int x, int y, int radius, int max_fcount, RGB RGB_table, linestringNode *tile_node, linestringNode *linestringRNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    double R0 = radius * Rz; 
    double R1 = R0 - sqrt(2) * Rz / 4.0;
    double R2 = R0 + sqrt(2) * Rz / 4.0;
    double Rz4 = 0.25 * Rz;
    int level = (radius == 0) ? 0 : int(radius + sqrt(2) / 4.0 - 0.5) + 1;
    // calcalate the BBox of tile <z,x,y>
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};

    int heat_matrix[256][256] = {0};

    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;

            // outer level
            if (i < level or i >= 256-level or j < level or j >= 256-level){
                for (int m = -level; m <= level; m++){
                    for (int n = -level; n <= level; n++){
                        double round_pix_x = pix_x + m * Rz;
                        double round_pix_y = pix_y + n * Rz;
                        if (sqrt(m * m + n * n) <= level){
                            linestringNode *round_pix_node = locPixNode(z, tile_box, round_pix_x, round_pix_y, tile_node, linestringRNode);
                            if (round_pix_node){
                                int round_fcount = round_pix_node->node_flen;
                                heat_matrix[i][j] += 255.0 * (1 - sqrt(m * m + n * n) / level) * (round_fcount - 0) / (max_fcount - 0);
                            }
                        }
                    }
                }
            }
            
            linestringNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
            if (pixel_node){
                int fcount = pixel_node->node_flen;

                for (int m = -level; m <= level; m++){
                    for (int n = -level; n <= level; n++){
                        if (sqrt(m * m + n * n) <= level){
                            // if (i+m < 256 and i+m >= 0 and j+n < 256 and j+n >= 0){
                                // heat_matrix[i+m][j+n] += 255.0 * (1 - sqrt(m * m + n * n) / level) * (fcount - 0) / (max_fcount - 0);
                            // }
                            if (i+m >= level and i+m < 256-level and j+n >= level and j+n < 256-level){
                                heat_matrix[i+m][j+n] += 255.0 * (1 - sqrt(m * m + n * n) / level) * (fcount - 0) / (max_fcount - 0);
                            }
                        }
                    }
                }
            }
        }
    }

    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            if (heat_matrix[i][j] > 0){
                int final_alpha = min(heat_matrix[i][j], 255);
                png_ptr[i][4 * j] =  RGB_table[final_alpha][0];     
                png_ptr[i][4 * j + 1] = RGB_table[final_alpha][1]; 
                png_ptr[i][4 * j + 2] = RGB_table[final_alpha][2]; 
                png_ptr[i][4 * j + 3] = final_alpha;
            }
            else{
                png_ptr[i][4 * j] =  0;    
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
        }
    }

}

// data-driven tile plotting
string data_tileRender_l_single(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, int R, int G, int B, float AD, float width, string data_string){
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt\n" + data_string;
    datasource_ptr line_ds = datasource_cache::instance().create(params);
    Map map(256, 256);
    map.set_srs(srs_merc);

    feature_type_style line_style;
    {
        rule r;
        {
            line_symbolizer line_sym;
            put(line_sym, keys::stroke, color(R, G, B));
            put(line_sym, keys::stroke_width, (width+1));
            put(line_sym, keys::stroke_linecap, ROUND_CAP);
            put(line_sym, keys::stroke_linejoin, ROUND_JOIN);
            put(line_sym, keys::stroke_opacity, AD);
            put(line_sym, keys::smooth, 1.0);
            r.append(std::move(line_sym));
        }
        line_style.add_rule(std::move(r));
    }
    map.insert_style("line_style", std::move(line_style));

    {
        layer lyr("line");
        lyr.set_datasource(line_ds);
        lyr.add_style("line_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    map.zoom_to_box(box2d<double>(tile_xMin, tile_yMin, tile_xMax, tile_yMax));
    image_rgba8 im(256, 256);
    agg_renderer<image_rgba8> ren(map, im);
    ren.apply();
    ostringstream im_os;
    save_to_stream(im, im_os, "png");
    return im_os.str();
}

string data_tileRender_l_range_single(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, double box[], int R, int G, int B, float AD, float width, string data_string){
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt\n" + data_string;
    datasource_ptr line_ds = datasource_cache::instance().create(params);
    Map map(256, 256);
    map.set_srs(srs_merc);

    feature_type_style line_style;
    {
        rule r;
        {
            line_symbolizer line_sym;
            put(line_sym, keys::stroke, color(R, G, B));
            put(line_sym, keys::stroke_width, (width+1));
            put(line_sym, keys::stroke_linecap, ROUND_CAP);
            put(line_sym, keys::stroke_linejoin, ROUND_JOIN);
            put(line_sym, keys::stroke_opacity, AD);
            put(line_sym, keys::smooth, 1.0);
            r.append(std::move(line_sym));
        }
        line_style.add_rule(std::move(r));
    }
    map.insert_style("line_style", std::move(line_style));

    {
        layer lyr("line");
        lyr.set_datasource(line_ds);
        lyr.add_style("line_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    map.zoom_to_box(box2d<double>(tile_xMin, tile_yMin, tile_xMax, tile_yMax));
    image_rgba8 im(256, 256);
    agg_renderer<image_rgba8> ren(map, im);
    ren.apply();

    double pix_Rz = (tile_xMax - tile_xMin) / 256;
    double box_xMin = box[0];
    double box_yMin = box[1];
    double box_xMax = box[2];
    double box_yMax = box[3];
    for (int w = 0; w < im.width(); w++){
        for (int h = 0; h < im.height(); h++){
            double pix_xMin = tile_xMin + h * pix_Rz;
            double pix_yMin = tile_yMax - (w + 1) * pix_Rz;

            if (!IfOverlapMBB(pix_xMin, pix_yMin, pix_Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                im(h, w) = color(255, 255, 255, 0).rgba();
            }
        }
    }

    ostringstream im_os;
    save_to_stream(im, im_os, "png");
    return im_os.str();
}

string data_tileRender_l_classify(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, string attri_name, vector<string> attri, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string data_string){
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt," + attri_name + "\n" + data_string;
    datasource_ptr line_ds = datasource_cache::instance().create(params);
    Map map(256, 256);
    map.set_srs(srs_merc);

    feature_type_style line_style;
    line_style.reserve(attri.size());
    for (int i = 0; i < attri.size(); i++){
        rule r;
        string style_describer = "[" + attri_name + "] = " + attri[i];
        r.set_filter(parse_expression(style_describer));
        {
            line_symbolizer line_sym;
            put(line_sym, keys::stroke, color(R[i], G[i], B[i]));
            put(line_sym, keys::stroke_width, (width[i]+1));
            put(line_sym, keys::stroke_linecap, ROUND_CAP);
            put(line_sym, keys::stroke_linejoin, ROUND_JOIN);
            put(line_sym, keys::stroke_opacity, AD[i]);
            put(line_sym, keys::smooth, 1.0);
            r.append(std::move(line_sym));
        }
        line_style.add_rule(std::move(r));
    }
    map.insert_style("line_style", std::move(line_style));

    {
        layer lyr("line");
        lyr.set_datasource(line_ds);
        lyr.add_style("line_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    map.zoom_to_box(box2d<double>(tile_xMin, tile_yMin, tile_xMax, tile_yMax));
    image_rgba8 im(256, 256);
    agg_renderer<image_rgba8> ren(map, im);
    ren.apply();
    ostringstream im_os;
    save_to_stream(im, im_os, "png");
    return im_os.str();
}

string data_tileRender_l_range_classify(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, double box[], string attri_name, vector<string> attri, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string data_string){
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt," + attri_name + "\n" + data_string;
    datasource_ptr line_ds = datasource_cache::instance().create(params);
    Map map(256, 256);
    map.set_srs(srs_merc);

    feature_type_style line_style;
    line_style.reserve(attri.size());
    for (int i = 0; i < attri.size(); i++){
        rule r;
        string style_describer = "[" + attri_name + "] = " + attri[i];
        r.set_filter(parse_expression(style_describer));
        {
            line_symbolizer line_sym;
            put(line_sym, keys::stroke, color(R[i], G[i], B[i]));
            put(line_sym, keys::stroke_width, (width[i]+1));
            put(line_sym, keys::stroke_linecap, ROUND_CAP);
            put(line_sym, keys::stroke_linejoin, ROUND_JOIN);
            put(line_sym, keys::stroke_opacity, AD[i]);
            put(line_sym, keys::smooth, 1.0);
            r.append(std::move(line_sym));
        }
        line_style.add_rule(std::move(r));
    }
    map.insert_style("line_style", std::move(line_style));

    {
        layer lyr("line");
        lyr.set_datasource(line_ds);
        lyr.add_style("line_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    map.zoom_to_box(box2d<double>(tile_xMin, tile_yMin, tile_xMax, tile_yMax));
    image_rgba8 im(256, 256);
    agg_renderer<image_rgba8> ren(map, im);
    ren.apply();

    double pix_Rz = (tile_xMax - tile_xMin) / 256;
    double box_xMin = box[0];
    double box_yMin = box[1];
    double box_xMax = box[2];
    double box_yMax = box[3];
    for (int w = 0; w < im.width(); w++){
        for (int h = 0; h < im.height(); h++){
            double pix_xMin = tile_xMin + h * pix_Rz;
            double pix_yMin = tile_yMax - (w + 1) * pix_Rz;

            if (!IfOverlapMBB(pix_xMin, pix_yMin, pix_Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                im(h, w) = color(255, 255, 255, 0).rgba();
            }
        }
    }

    ostringstream im_os;
    save_to_stream(im, im_os, "png");
    return im_os.str();
}

string data_tileRender_l_rule(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, string attri_name, vector<string> filters, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string data_string){
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt," + attri_name + "\n" + data_string;
    datasource_ptr line_ds = datasource_cache::instance().create(params);
    Map map(256, 256);
    map.set_srs(srs_merc);

    feature_type_style line_style;
    line_style.reserve(filters.size());
    for (int i = 0; i < filters.size(); i++){
        rule r;
        int ix;
        string filter_i = filters[i];

        if ((ix = filter_i.find("%3E=")) != string::npos){
            filter_i.replace(0, 4, ">=");
        }
        else if ((ix = filter_i.find("%3C=")) != string::npos){
            filter_i.replace(0, 4, "<=");
        }
        else if ((ix = filter_i.find("%3E")) != string::npos){
            filter_i.replace(0, 3, ">");
        }
        else if ((ix = filter_i.find("%3C")) != string::npos){
            filter_i.replace(0, 3, "<");
        }
        string style_describer = "[" + attri_name + "] " + filter_i;


        r.set_filter(parse_expression(style_describer));
        {
            line_symbolizer line_sym;
            put(line_sym, keys::stroke, color(R[i], G[i], B[i]));
            put(line_sym, keys::stroke_width, (width[i]+1));
            put(line_sym, keys::stroke_linecap, ROUND_CAP);
            put(line_sym, keys::stroke_linejoin, ROUND_JOIN);
            put(line_sym, keys::stroke_opacity, AD[i]);
            put(line_sym, keys::smooth, 1.0);
            r.append(std::move(line_sym));
        }
        line_style.add_rule(std::move(r));
    }
    map.insert_style("line_style", std::move(line_style));

    {
        layer lyr("line");
        lyr.set_datasource(line_ds);
        lyr.add_style("line_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    map.zoom_to_box(box2d<double>(tile_xMin, tile_yMin, tile_xMax, tile_yMax));
    image_rgba8 im(256, 256);
    agg_renderer<image_rgba8> ren(map, im);
    ren.apply();
    ostringstream im_os;
    save_to_stream(im, im_os, "png");
    return im_os.str();
}

string data_tileRender_l_range_rule(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, double box[], string attri_name, vector<string> filters, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string data_string){
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt," + attri_name + "\n" + data_string;
    datasource_ptr line_ds = datasource_cache::instance().create(params);
    Map map(256, 256);
    map.set_srs(srs_merc);

    feature_type_style line_style;
    line_style.reserve(filters.size());
    for (int i = 0; i < filters.size(); i++){
        rule r;
        int ix;
        string filter_i = filters[i];

        if ((ix = filter_i.find("%3E=")) != string::npos){
            filter_i.replace(0, 4, ">=");
        }
        else if ((ix = filter_i.find("%3C=")) != string::npos){
            filter_i.replace(0, 4, "<=");
        }
        else if ((ix = filter_i.find("%3E")) != string::npos){
            filter_i.replace(0, 3, ">");
        }
        else if ((ix = filter_i.find("%3C")) != string::npos){
            filter_i.replace(0, 3, "<");
        }
        string style_describer = "[" + attri_name + "] " + filter_i;


        r.set_filter(parse_expression(style_describer));
        {
            line_symbolizer line_sym;
            put(line_sym, keys::stroke, color(R[i], G[i], B[i]));
            put(line_sym, keys::stroke_width, (width[i]+1));
            put(line_sym, keys::stroke_linecap, ROUND_CAP);
            put(line_sym, keys::stroke_linejoin, ROUND_JOIN);
            put(line_sym, keys::stroke_opacity, AD[i]);
            put(line_sym, keys::smooth, 1.0);
            r.append(std::move(line_sym));
        }
        line_style.add_rule(std::move(r));
    }
    map.insert_style("line_style", std::move(line_style));

    {
        layer lyr("line");
        lyr.set_datasource(line_ds);
        lyr.add_style("line_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    map.zoom_to_box(box2d<double>(tile_xMin, tile_yMin, tile_xMax, tile_yMax));
    image_rgba8 im(256, 256);
    agg_renderer<image_rgba8> ren(map, im);
    ren.apply();

    double pix_Rz = (tile_xMax - tile_xMin) / 256;
    double box_xMin = box[0];
    double box_yMin = box[1];
    double box_xMax = box[2];
    double box_yMax = box[3];
    for (int w = 0; w < im.width(); w++){
        for (int h = 0; h < im.height(); h++){
            double pix_xMin = tile_xMin + h * pix_Rz;
            double pix_yMin = tile_yMax - (w + 1) * pix_Rz;

            if (!IfOverlapMBB(pix_xMin, pix_yMin, pix_Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                im(h, w) = color(255, 255, 255, 0).rgba();
            }
        }
    }

    ostringstream im_os;
    save_to_stream(im, im_os, "png");
    return im_os.str();
}

// diaplay-driven visualization
string disDrivenPlot(int z, int x, int y, string style_describer, linestringNode *RNode){
    png_bytep * tile_ptr = (png_bytep *) malloc(256 * sizeof(png_bytep));
    for (int n = 0; n < TILE_SIZE; n++)
        tile_ptr[n] = (png_bytep) malloc(1024);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr); 
    FILE *temp_png = tmpfile();
    png_init_io(png_ptr, temp_png); 
    png_set_IHDR(png_ptr, info_ptr, TILE_SIZE, TILE_SIZE, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    linestringNode *tile_node = locTileNode(z, x, y, RNode);
    if (tile_node != nullptr){
        vector<string> attri_list;
        boost::split(attri_list, style_describer, boost::is_any_of( "-" ), boost::token_compress_on);
        string style_type = attri_list[0];
        // single visualization
        if (style_type == "single"){
            vector<string> style_list;
            boost::split(style_list, attri_list[1], boost::is_any_of( "," ), boost::token_compress_on);
            int R = stoi(style_list[0]);
            int G = stoi(style_list[1]);
            int B = stoi(style_list[2]);
            float AD = stof(style_list[3]);
            float width = stof(style_list[4]);
            dis_tileRender_single(z, x, y, R, G, B, AD, width, tile_node, RNode, tile_ptr);
        }     
        // classified visualization   
        else if (style_type == "classified"){
            vector<string> attri;
            vector<int> R, G, B;
            vector<float> AD, width;
            string attri_name = attri_list[1];
            string attri_type = attri_list[2];
            for (int i = 3; i < attri_list.size(); i++){
                vector<string> style_list;
                boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                attri.push_back(style_list[0]);
                R.push_back(stoi(style_list[1]));
                G.push_back(stoi(style_list[2]));
                B.push_back(stoi(style_list[3]));
                AD.push_back(stof(style_list[4]));
                width.push_back(stof(style_list[5]));
            }
            dis_tileRender_classify(z, x, y, attri, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
        }
        // ruled visualization
        else if (style_type == "ruled"){
            vector<string> filters;
            vector<int> R, G, B;
            vector<float> AD, width;
            string attri_name = attri_list[1];
            string attri_type = attri_list[2];
            for (int i = 3; i < attri_list.size(); i++){
                vector<string> style_list;
                boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                filters.push_back(style_list[0]);
                R.push_back(stoi(style_list[1]));
                G.push_back(stoi(style_list[2]));
                B.push_back(stoi(style_list[3]));
                AD.push_back(stof(style_list[4]));
                width.push_back(stof(style_list[5]));
            }
            dis_tileRender_rule(z, x, y, filters, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
        }
        // heatmap visualization
        else if (style_type == "heat"){
            string gradient_type = attri_list[1];
            int radius = stoi(attri_list[2]);
            int max_fcount = 0;
            RNode->getMaxFcount(z, max_fcount);
            RGB RGB_table = initRGBTable(gradient_type);
            dis_tileRender_heat(z, x, y, radius, max_fcount, RGB_table, tile_node, RNode, tile_ptr);
        }
        else{
            # pragma omp parallel for num_threads(4)
            for (int i = 0; i < 256; i++){
                for (int j = 0; j < 256 * 4; j++){
                    tile_ptr[i][j] = 0; 
                }
            }
        }
    }
    else{
        # pragma omp parallel for num_threads(4)
        for (int i = 0; i < 256; i++){
            for (int j = 0; j < 256 * 4; j++){
                tile_ptr[i][j] = 0; 
            }
        }
    }

    png_write_image(png_ptr, tile_ptr);   
    png_write_end(png_ptr, NULL);
    for (int k = 0; k < 256; k++)
        free(tile_ptr[k]);
    free(tile_ptr);
    vector<char> pos;
    fseek(temp_png, 0, SEEK_END);
    long size = ftell(temp_png);
    rewind(temp_png);
    pos.resize(size);
    fread(&pos[0], 1, size, temp_png);
    fclose(temp_png);
    string pos_tostr = string(pos.begin(), pos.end());
    return pos_tostr;
}

string disDrivenPlot_range(int z, int x, int y, double box[], string style_describer, linestringNode *RNode){
    png_bytep * tile_ptr = (png_bytep *) malloc(256 * sizeof(png_bytep));
    for (int n = 0; n < TILE_SIZE; n++)
        tile_ptr[n] = (png_bytep) malloc(1024);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr); 
    FILE *temp_png = tmpfile();
    png_init_io(png_ptr, temp_png); 
    png_set_IHDR(png_ptr, info_ptr, TILE_SIZE, TILE_SIZE, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    linestringNode *tile_node = locTileNode(z, x, y, RNode);
    if (tile_node != nullptr){
        double tile_Rz = L / pow(2, z-1);
        double tile_xMin = x * tile_Rz - L;
        double tile_yMin = L - (y + 1) * tile_Rz;
        double box_xMin = box[0];
        double box_yMin = box[1];
        double box_xMax = box[2];
        double box_yMax = box[3];
        vector<string> attri_list;
        boost::split(attri_list, style_describer, boost::is_any_of( "-" ), boost::token_compress_on);
        string style_type = attri_list[0];
        // single visualization
        if (style_type == "single"){
            vector<string> style_list;
            boost::split(style_list, attri_list[1], boost::is_any_of( "," ), boost::token_compress_on);
            int R = stoi(style_list[0]);
            int G = stoi(style_list[1]);
            int B = stoi(style_list[2]);
            float AD = stof(style_list[3]);
            float width = stof(style_list[4]);
            if (IfContain(box_xMin, box_yMin, box_xMax, box_yMax, tile_xMin, tile_yMin, tile_xMin+tile_Rz, tile_yMin+tile_Rz)){
                dis_tileRender_single(z, x, y, R, G, B, AD, width, tile_node, RNode, tile_ptr);
            }
            else{
                dis_tileRender_range_single(z, x, y, box, R, G, B, AD, width, tile_node, RNode, tile_ptr);
            }
            
        }     
        else{
            vector<string> filters;
            vector<int> R, G, B;
            vector<float> AD, width;
            string attri_name = attri_list[1];
            string attri_type = attri_list[2];
            for (int i = 3; i < attri_list.size(); i++){
                vector<string> style_list;
                boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                filters.push_back(style_list[0]);
                R.push_back(stoi(style_list[1]));
                G.push_back(stoi(style_list[2]));
                B.push_back(stoi(style_list[3]));
                AD.push_back(stof(style_list[4]));
                width.push_back(stof(style_list[5]));
            }
            // classified visualization   
            if (style_type == "classified"){
                if (IfContain(box_xMin, box_yMin, box_xMax, box_yMax, tile_xMin, tile_yMin, tile_xMin+tile_Rz, tile_yMin+tile_Rz)){
                    dis_tileRender_classify(z, x, y, filters, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
                }
                else{
                    dis_tileRender_range_classify(z, x, y, box, filters, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
                }
            }
            
            // ruled visualization
            else if (style_type == "ruled"){
                if (IfContain(box_xMin, box_yMin, box_xMax, box_yMax, tile_xMin, tile_yMin, tile_xMin+tile_Rz, tile_yMin+tile_Rz)){
                    dis_tileRender_rule(z, x, y, filters, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
                }
                else{
                    dis_tileRender_range_rule(z, x, y, box, filters, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
                }
            }
        }
    }
    else{
        # pragma omp parallel for num_threads(4)
        for (int i = 0; i < 256; i++){
            for (int j = 0; j < 256 * 4; j++){
                tile_ptr[i][j] = 0; 
            }
        }
    }

    png_write_image(png_ptr, tile_ptr);   
    png_write_end(png_ptr, NULL);
    for (int k = 0; k < 256; k++)
        free(tile_ptr[k]);
    free(tile_ptr);
    vector<char> pos;
    fseek(temp_png, 0, SEEK_END);
    long size = ftell(temp_png);
    rewind(temp_png);
    pos.resize(size);
    fread(&pos[0], 1, size, temp_png);
    fclose(temp_png);
    string pos_tostr = string(pos.begin(), pos.end());
    return pos_tostr;
}

// data-driven visualization
string dataDrivenPlot(int z, int x, int y, string style_describer, linestringNode *RNode, short switch_level){
    string img_str;
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    double tile_xCenter = tile_xMin + (tile_xMax - tile_xMin) / 2;
    double tile_yCenter = tile_yMin + (tile_yMax - tile_yMin) / 2;
    // determine whether the tile needs to be drawn
    short linestring_treeLevel = switch_level + 9;
    linestringNode *rtree_node = locRtreeNode(z, tile_xCenter, tile_yCenter, RNode, switch_level);
    // non-blank tile
    if (rtree_node != nullptr){
        // get the features by spatial range retrieve
        vector<LINESTRING_ITEM> objs;
        Box box(POINT(tile_xMin, tile_yMin), POINT(tile_xMax, tile_yMax));
        rtree_node->Rtree.query(bgi::intersects(box), back_inserter(objs));
        string data_string = "";
        // plot features
        vector<string> attri_list;
        boost::split(attri_list, style_describer, boost::is_any_of( "-" ), boost::token_compress_on);
        string style_type = attri_list[0];
        // single visualization
        if (style_type == "single"){
            vector<string> style_list;
            boost::split(style_list, attri_list[1], boost::is_any_of( "," ), boost::token_compress_on);
            int R = stoi(style_list[0]);
            int G = stoi(style_list[1]);
            int B = stoi(style_list[2]);
            float AD = stof(style_list[3]);
            float width = stof(style_list[4]);
            for (int i = 0; i < objs.size(); i++){
                LINE line_coord = get<1>(objs[i]);
                data_string += "\"LINESTRING(";
                for (int j = 0; j < line_coord.size()-2; j+=2){
                    data_string += (to_string(line_coord[j]) + " " + to_string(line_coord[j+1]) + ",");
                }
                data_string += (to_string(line_coord[line_coord.size()-2]) + " " + to_string(line_coord.back()) + ")\"\n");

                // data_string += ("\"" + get<1>(objs[i]) + "\"\n");
            }
            img_str = data_tileRender_l_single(tile_xMin, tile_yMin, tile_xMax, tile_yMax, R, G, B, AD, width, data_string);
        }
        else{
            vector<int> R, G, B;
            vector<float> AD, width;
            string attri_name = attri_list[1];
            string attri_type = attri_list[2];
            for (int i = 0; i < objs.size(); i++){
                LINE line_coord = get<1>(objs[i]);
                linestringAttri attri = get<2>(objs[i]);
                data_string += "\"LINESTRING(";
                for (int j = 0; j < line_coord.size()-2; j+=2){
                    data_string += (to_string(line_coord[j]) + " " + to_string(line_coord[j+1]) + ",");
                }
                data_string += (to_string(line_coord[line_coord.size()-2]) + " " + to_string(line_coord.back()) + ")\",");
                data_string += (attri.field_string[attri_name] + "\n");

                // data_string += ("\"" + get<1>(objs[i]) + "\"," + attri.field[attri_name] + "\n");
            }

            // classified visualization   
            if (style_type == "classified"){
                vector<string> attri;
                for (int i = 3; i < attri_list.size(); i++){
                    vector<string> style_list;
                    boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                    attri.push_back(style_list[0]);
                    R.push_back(stoi(style_list[1]));
                    G.push_back(stoi(style_list[2]));
                    B.push_back(stoi(style_list[3]));
                    AD.push_back(stof(style_list[4]));
                    width.push_back(stof(style_list[5]));
                }
                img_str = data_tileRender_l_classify(tile_xMin, tile_yMin, tile_xMax, tile_yMax, attri_name, attri, R, G, B, AD, width, data_string);
            }
            // ruled visualization
            else if (style_type == "ruled"){
                vector<string> rules;
                for (int i = 3; i < attri_list.size(); i++){
                    vector<string> style_list;
                    boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                    rules.push_back(style_list[0]);
                    R.push_back(stoi(style_list[1]));
                    G.push_back(stoi(style_list[2]));
                    B.push_back(stoi(style_list[3]));
                    AD.push_back(stof(style_list[4]));
                    width.push_back(stof(style_list[5]));
                }
                img_str = data_tileRender_l_rule(tile_xMin, tile_yMin, tile_xMax, tile_yMax, attri_name, rules, R, G, B, AD, width, data_string);
            }
            else{
                Map map(256, 256);
                const color back_color(255, 255, 255, 0);
                map.set_background(back_color);
                image_rgba8 im(256, 256);
                agg_renderer<image_rgba8> ren(map, im);
                ren.apply();
                ostringstream im_os;
                save_to_stream(im, im_os, "png");
                img_str = im_os.str();
            }
        }
    }
    // blank tile
    else{
        Map map(256, 256);
        const color back_color(255, 255, 255, 0);
        map.set_background(back_color);
        image_rgba8 im(256, 256);
        agg_renderer<image_rgba8> ren(map, im);
        ren.apply();
        ostringstream im_os;
        save_to_stream(im, im_os, "png");
        img_str = im_os.str();
    }
    return img_str;
}

string dataDrivenPlot_range(int z, int x, int y, double range_box[], string style_describer, linestringNode *RNode, short switch_level){
    string img_str;
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    double tile_xCenter = tile_xMin + (tile_xMax - tile_xMin) / 2;
    double tile_yCenter = tile_yMin + (tile_yMax - tile_yMin) / 2;
    // determine whether the tile needs to be drawn
    short linestring_treeLevel = switch_level + 9;
    linestringNode *rtree_node = locRtreeNode(z, tile_xCenter, tile_yCenter, RNode, switch_level);
    // non-blank tile
    if (rtree_node != nullptr){
        double box_xMin = range_box[0];
        double box_yMin = range_box[1];
        double box_xMax = range_box[2];
        double box_yMax = range_box[3];
        // get the features by spatial range retrieve
        vector<LINESTRING_ITEM> objs;
        Box box(POINT(tile_xMin, tile_yMin), POINT(tile_xMax, tile_yMax));
        rtree_node->Rtree.query(bgi::intersects(box), back_inserter(objs));
        string data_string = "";
        // plot features
        vector<string> attri_list;
        boost::split(attri_list, style_describer, boost::is_any_of( "-" ), boost::token_compress_on);
        string style_type = attri_list[0];
        // single visualization
        if (style_type == "single"){
            vector<string> style_list;
            boost::split(style_list, attri_list[1], boost::is_any_of( "," ), boost::token_compress_on);
            int R = stoi(style_list[0]);
            int G = stoi(style_list[1]);
            int B = stoi(style_list[2]);
            float AD = stof(style_list[3]);
            float width = stof(style_list[4]);
            for (int i = 0; i < objs.size(); i++){
                LINE line_coord = get<1>(objs[i]);
                data_string += "\"LINESTRING(";
                for (int j = 0; j < line_coord.size()-2; j+=2){
                    data_string += (to_string(line_coord[j]) + " " + to_string(line_coord[j+1]) + ",");
                }
                data_string += (to_string(line_coord[line_coord.size()-2]) + " " + to_string(line_coord.back()) + ")\"\n");
            }

            if (IfContain(box_xMin, box_yMin, box_xMax, box_yMax, tile_xMin, tile_yMin, tile_xMax, tile_yMax)){
                img_str = data_tileRender_l_single(tile_xMin, tile_yMin, tile_xMax, tile_yMax, R, G, B, AD, width, data_string);
            }
            else{
                img_str = data_tileRender_l_range_single(tile_xMin, tile_yMin, tile_xMax, tile_yMax, range_box, R, G, B, AD, width, data_string);
            }
        }
        else{
            vector<int> R, G, B;
            vector<float> AD, width;
            string attri_name = attri_list[1];
            string attri_type = attri_list[2];
            for (int i = 0; i < objs.size(); i++){
                LINE line_coord = get<1>(objs[i]);
                linestringAttri attri = get<2>(objs[i]);
                data_string += "\"LINESTRING(";
                for (int j = 0; j < line_coord.size()-2; j+=2){
                    data_string += (to_string(line_coord[j]) + " " + to_string(line_coord[j+1]) + ",");
                }
                data_string += (to_string(line_coord[line_coord.size()-2]) + " " + to_string(line_coord.back()) + ")\",");
                data_string += (attri.field_string[attri_name] + "\n");
            }

            // classified visualization   
            if (style_type == "classified"){
                vector<string> attri;
                for (int i = 3; i < attri_list.size(); i++){
                    vector<string> style_list;
                    boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                    attri.push_back(style_list[0]);
                    R.push_back(stoi(style_list[1]));
                    G.push_back(stoi(style_list[2]));
                    B.push_back(stoi(style_list[3]));
                    AD.push_back(stof(style_list[4]));
                    width.push_back(stof(style_list[5]));
                }

                if (IfContain(box_xMin, box_yMin, box_xMax, box_yMax, tile_xMin, tile_yMin, tile_xMax, tile_yMax)){
                    img_str = data_tileRender_l_classify(tile_xMin, tile_yMin, tile_xMax, tile_yMax, attri_name, attri, R, G, B, AD, width, data_string);
                }
                else{
                    img_str = data_tileRender_l_range_classify(tile_xMin, tile_yMin, tile_xMax, tile_yMax, range_box, attri_name, attri, R, G, B, AD, width, data_string);
                }
            }
            // ruled visualization
            else if (style_type == "ruled"){
                vector<string> rules;
                for (int i = 3; i < attri_list.size(); i++){
                    vector<string> style_list;
                    boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                    rules.push_back(style_list[0]);
                    R.push_back(stoi(style_list[1]));
                    G.push_back(stoi(style_list[2]));
                    B.push_back(stoi(style_list[3]));
                    AD.push_back(stof(style_list[4]));
                    width.push_back(stof(style_list[5]));
                }

                if (IfContain(box_xMin, box_yMin, box_xMax, box_yMax, tile_xMin, tile_yMin, tile_xMax, tile_yMax)){
                    img_str = data_tileRender_l_rule(tile_xMin, tile_yMin, tile_xMax, tile_yMax, attri_name, rules, R, G, B, AD, width, data_string);
                }
                else{
                    img_str = data_tileRender_l_range_rule(tile_xMin, tile_yMin, tile_xMax, tile_yMax, range_box, attri_name, rules, R, G, B, AD, width, data_string);
                }
            }
            else{
                Map map(256, 256);
                const color back_color(255, 255, 255, 0);
                map.set_background(back_color);
                image_rgba8 im(256, 256);
                agg_renderer<image_rgba8> ren(map, im);
                ren.apply();
                ostringstream im_os;
                save_to_stream(im, im_os, "png");
                img_str = im_os.str();
            }
        }
    }
    // blank tile
    else{
        Map map(256, 256);
        const color back_color(255, 255, 255, 0);
        map.set_background(back_color);
        image_rgba8 im(256, 256);
        agg_renderer<image_rgba8> ren(map, im);
        ren.apply();
        ostringstream im_os;
        save_to_stream(im, im_os, "png");
        img_str = im_os.str();
    }
    return img_str;
}



string fieldStatistic(double box_xMin, double box_yMin, double box_xMax, double box_yMax, string style_describer, short switch_level, linestringNode *RNode){
    string field_statistic = "";
    map<string, double> count_list;
    vector<string> attri_list;
    vector<string> style_list;
    boost::split(attri_list, style_describer, boost::is_any_of( "-" ), boost::token_compress_on);
    string style_type = attri_list[0];

    if (style_type == "single"){
        linestringNode *min_WNode = RNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        
        if (min_WNode){
            min_WNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_list);
            boost::split(style_list, attri_list[1], boost::is_any_of( "," ), boost::token_compress_on);
            if (count_list.size()){
                field_statistic = "Count:" + to_string(count_list["length"]);
            }
            else{
                field_statistic = "Count:0";
            }
        }
        else{
            field_statistic = "Count:0";
        }
        
    }
    if (style_type == "classified"){
        string attri_name = attri_list[1];
        linestringNode *min_WNode = RNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        if (min_WNode){
            min_WNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, attri_name, count_list);
            for (int i = 3; i < attri_list.size(); i++){
                boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                string attri = style_list[0];
                double field_value = count_list[attri];
                field_statistic += (attri + ":" + to_string(field_value) + ",");
            }
            field_statistic.pop_back();
        }
        
    }
    else if (style_type == "ruled"){
        string attri_name = attri_list[1];
        linestringNode *min_WNode = RNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        if (min_WNode){
            vector<string> rules, R, G, B, AD, width;
            for (int i = 3; i < attri_list.size(); i++){
                boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                rules.push_back(style_list[0]);
            }
            min_WNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, attri_name, rules, count_list);
            for (int j = 0; j < rules.size(); j++){
                double field_value = count_list[rules[j]];
                field_statistic += (rules[j] + ":" + to_string(field_value) + ",");
            }
            field_statistic.pop_back();
        }
        
    }

    return field_statistic;
}