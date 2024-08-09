#include <iostream>
#include <fstream>
#include "png.h"
#include "polygonPlot.hpp"
#include "treeNode.hpp"
#include "pointPlot.hpp"
#include "linestringPlot.hpp"
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



// determine whether the node correspongding to the pixel exists
char nodeIfExist(int z, double tile_box[], double pix_x, double pix_y, polygonNode *seek){
    double min_x = tile_box[0];
    double min_y = tile_box[1];
    double max_x = tile_box[2];
    double max_y = tile_box[3];
    double mid_x, mid_y;
    char node_topo = 'n';

    int divided_num = 0;
    int did_num = 8 + z - seek->node_level;
    while (divided_num < did_num){
        if (seek->topo_type == 'w'){
            // node_topo = 'w';
            return 'w';
        }
        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // pixel in LU sub-region
        if (pix_x < mid_x && pix_y >= mid_y){
            if (seek->LUNode == nullptr){
                // if (node_topo == 'w')
                //     return 'w';
                return 'n';
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // pixel in RU sub-region
        else if (pix_x >= mid_x && pix_y >= mid_y){
            if (seek->RUNode == nullptr){
                // if (node_topo == 'w')
                //     return 'w';
                return 'n';
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // pixel in LB sub-region
        else if (pix_x < mid_x && pix_y < mid_y){
            if (seek->LBNode == nullptr){
                // if (node_topo == 'w')
                //     return 'w';
                return 'n';
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // pixel in RB sub-region
        else if (pix_x >= mid_x && pix_y < mid_y){
            if (seek->RBNode == nullptr){
                // if (node_topo == 'w')
                //     return 'w';
                return 'n';
            }
            seek = seek->RBNode;
            min_x = mid_x;
            max_y = mid_y;
        }
        divided_num++;
    }    

    return seek->topo_type;
}


polygonNode* nodeIfExist_attri(int z, double tile_box[], double pix_x, double pix_y, polygonNode *seek){
    float min_x = tile_box[0];
    float min_y = tile_box[1];
    float max_x = tile_box[2];
    float max_y = tile_box[3];
    double mid_x, mid_y;

    int divided_num = 0;

    int did_num = 8;
    while (divided_num < did_num){
        if (seek->topo_type == 'w'){
            return seek;
        }
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
        }
        divided_num++;
    }

    return seek;
}

// determine whether the node correspongding to the tile<z,x,y> exists
polygonNode * locTileNode(int z, int x, int y, polygonNode *seek){
    float tile_Rz = L / pow(2, z-1);
    float tile_x = (x + 0.5) * tile_Rz - L;
    float tile_y = L - (y + 0.5) * tile_Rz;
    float min_x = -L;
    float max_x = L;
    float min_y = -L;
    float max_y = L;
    float mid_x, mid_y;
    int divided_num = 0;

    while (divided_num < z){
        if (seek->topo_type == 'w'){
            return seek;
        }

        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // tile_center_pt in LU sub-region
        if (tile_x < mid_x && tile_y >= mid_y){
            if (seek->LUNode == nullptr){
                return nullptr;
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in RU sub-region
        else if (tile_x >= mid_x && tile_y >= mid_y){
            if (seek->RUNode == nullptr){
                return nullptr;
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in LB sub-region
        else if (tile_x < mid_x && tile_y < mid_y){
            if (seek->LBNode == nullptr){
                return nullptr;
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // tile_center_pt in RB sub-region
        else if (tile_x >= mid_x && tile_y < mid_y){
            if (seek->RBNode == nullptr){
                return nullptr;
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
pair<polygonNode*, polygonNode*> locRtreeNode(int did_num, int tile_x, int tile_y, polygonNode *seek, short switch_level){
    float min_x = -L;
    float max_x = L;
    float min_y = -L;
    float max_y = L;
    float mid_x, mid_y;
    int divided_num = 0;
    polygonNode *middle_node = nullptr;

    if (did_num > switch_level + 9)
        did_num = switch_level + 9;

    while (divided_num < did_num){
        if (seek->topo_type == 'w'){
            return make_pair(seek, middle_node);
        }

        mid_x = (min_x + max_x) / 2;
        mid_y = (min_y + max_y) / 2;
        // tile_center_pt in LU sub-region
        if (tile_x < mid_x && tile_y >= mid_y){
            if (seek->LUNode == nullptr){
                return make_pair(nullptr, nullptr);
            }
            seek = seek->LUNode;
            max_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in RU sub-region
        else if (tile_x >= mid_x && tile_y >= mid_y){
            if (seek->RUNode == nullptr){
                return make_pair(nullptr, nullptr);
            } 
            seek = seek->RUNode;
            min_x = mid_x;
            min_y = mid_y;
        }
        // tile_center_pt in LB sub-region
        else if (tile_x < mid_x && tile_y < mid_y){
            if (seek->LBNode == nullptr){
                return make_pair(nullptr, nullptr);
            }
            seek = seek->LBNode;
            max_x = mid_x;
            max_y = mid_y;
        }
        // tile_center_pt in RB sub-region
        else if (tile_x >= mid_x && tile_y < mid_y){
            if (seek->RBNode == nullptr){
                return make_pair(nullptr, nullptr);
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

    return make_pair(seek, middle_node);
}


// determine whether the node correspongding to the nearest pixel exists
polygonNode * locPixNode(int z, double tile_box[], double pix_x, double pix_y, polygonNode *seek, polygonNode *rootNode){
    double min_x = tile_box[0];
    double min_y = tile_box[1];
    double max_x = tile_box[2];
    double max_y = tile_box[3];
    double mid_x, mid_y;

    if (abs(pix_x) <= L && abs(pix_y) <= L){
        if (pix_x < min_x || pix_x > max_x || pix_y < min_y || pix_y > max_y){
            seek = rootNode;
            min_x = -L;
            max_x = L;
            min_y = -L;
            max_y = L; 
        }
    }
    else{
        return nullptr;
    }

    int level = 8 + z - seek->node_level;
    int divided_num = 0;

    while (divided_num < level){
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
        }
        divided_num++;
    }

    if (seek->topo_type == 'w')
        return nullptr;

    return seek;
}



// get the child node of the node which is nearest to the pixel 
double *getNearestCNode(polygonNode *seek, double pix_x, double pix_y, double near_x, double near_y, double Rz){
    double dist = 0;
    double node_x, node_y;
    if (seek->LUNode != nullptr && seek->LUNode->topo_type == 'i'){
        double node_x1 = near_x - Rz;
        double node_y1 = near_y + Rz;
        double dist1 = sqrt(pow(node_x1 - pix_x, 2) + pow(node_y1 - pix_y, 2));
        if (dist == 0){
            dist = dist1;
            node_x = node_x1;
            node_y = node_y1;
        }
    }
    if (seek->RUNode != nullptr && seek->RUNode->topo_type == 'i'){
        double node_x2 = near_x + Rz;
        double node_y2 = near_y + Rz;
        double dist2 = sqrt(pow(node_x2 - pix_x, 2) + pow(node_y2 - pix_y, 2));
        if (dist == 0 || dist2 < dist){
            dist = dist2;
            node_x = node_x2;
            node_y = node_y2;
        }
    }
    if (seek->LBNode != nullptr && seek->LBNode->topo_type == 'i'){
        double node_x3 = near_x - Rz;
        double node_y3 = near_y - Rz;
        double dist3 = sqrt(pow(node_x3 - pix_x, 2) + pow(node_y3 - pix_y, 2));
        if (dist == 0 || dist3 < dist){
            dist = dist3;
            node_x = node_x3;
            node_y = node_y3;
        }
    }
    if (seek->RBNode != nullptr && seek->RBNode->topo_type == 'i'){
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
void dis_tileRender_single(int z, int x, int y, int R, int G, int B, float AD, float width, int f_R, int f_G, int f_B, float f_AD, polygonNode *tile_node, polygonNode *RNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    double R0 = width * Rz; 
    double R1 = R0 - sqrt(2) * Rz / 4.0;
    double R2 = R0 + sqrt(2) * Rz / 4.0;
    double Rz4 = 0.25 * Rz;
    int level = (width == 0) ? 0 : ceil(width + sqrt(2) / 4.0 - 0.5);
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
            char node_topoType = nodeIfExist(z, tile_box, pix_x, pix_y, tile_node);
            if (node_topoType == 'i'){
                png_ptr[i][4 * j] = R;
                png_ptr[i][4 * j + 1] = G;
                png_ptr[i][4 * j + 2] = B;
                png_ptr[i][4 * j + 3] = 255 * AD;
                continue;
            }
            else if (node_topoType == 'w'){
                png_ptr[i][4 * j] = f_R;
                png_ptr[i][4 * j + 1] = f_G;
                png_ptr[i][4 * j + 2] = f_B;
                png_ptr[i][4 * j + 3] = 255 * f_AD;
            }
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
            
            short final_value = 0;
            short pix_value = 0;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (abs(m) == l || abs(n) == l){
                            double nearPix_x = pix_x + m * Rz;
                            double nearPix_y = pix_y + n * Rz;
                            polygonNode *node_floor = locPixNode(z, tile_box, nearPix_x, nearPix_y, tile_node, RNode);
                            if (node_floor != nullptr){
                                    double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, nearPix_x, nearPix_y, Rz / 4);
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
            if (final_value){
                png_ptr[i][4 * j] = R;
                png_ptr[i][4 * j + 1] = G;
                png_ptr[i][4 * j + 2] = B;
                png_ptr[i][4 * j + 3] = AD * final_value * 64 - 1;
            }
            if ((final_value == 1) and node_topoType == 'w'){
                png_ptr[i][4 * j] = f_R;
                png_ptr[i][4 * j + 1] = f_G;
                png_ptr[i][4 * j + 2] = f_B;
                png_ptr[i][4 * j + 3] = f_AD * final_value * 64 - 1;
            }
        }
    }
}

void dis_tileRender_range_single(int z, int x, int y, double box[], int R, int G, int B, float AD, float width, int f_R, int f_G, int f_B, float f_AD, polygonNode *tile_node, polygonNode *RNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    double R0 = width * Rz; 
    double R1 = R0 - sqrt(2) * Rz / 4.0;
    double R2 = R0 + sqrt(2) * Rz / 4.0;
    double Rz4 = 0.25 * Rz;
    int level = (width == 0) ? 0 : ceil(width + sqrt(2) / 4.0 - 0.5);
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
            char node_topoType = nodeIfExist(z, tile_box, pix_x, pix_y, tile_node);
            if (IfOverlapMBB(pix_xMin, pix_yMin, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                // determine whether the node correspongding to the pixel exists
                if (node_topoType == 'i'){
                    png_ptr[i][4 * j] = R;
                    png_ptr[i][4 * j + 1] = G;
                    png_ptr[i][4 * j + 2] = B;
                    png_ptr[i][4 * j + 3] = 255 * AD;
                    continue;
                }
                else if (node_topoType == 'w'){
                    png_ptr[i][4 * j] = f_R;
                    png_ptr[i][4 * j + 1] = f_G;
                    png_ptr[i][4 * j + 2] = f_B;
                    png_ptr[i][4 * j + 3] = 255 * f_AD;
                }
                else{
                    png_ptr[i][4 * j] = 0;
                    png_ptr[i][4 * j + 1] = 0;
                    png_ptr[i][4 * j + 2] = 0;
                    png_ptr[i][4 * j + 3] = 0;
                }
            }
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }
            
            short final_value = 0;
            short pix_value = 0;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (abs(m) == l || abs(n) == l){
                            if (IfOverlapMBB(pix_xMin + m * Rz, pix_yMin + n * Rz, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                                double nearPix_x = pix_x + m * Rz;
                                double nearPix_y = pix_y + n * Rz;
                                polygonNode *node_floor = locPixNode(z, tile_box, nearPix_x, nearPix_y, tile_node, RNode);
                                if (node_floor != nullptr){
                                        double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, nearPix_x, nearPix_y, Rz / 4);
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
            if (final_value){
                png_ptr[i][4 * j] = R;
                png_ptr[i][4 * j + 1] = G;
                png_ptr[i][4 * j + 2] = B;
                png_ptr[i][4 * j + 3] = AD * final_value * 64 - 1;
            }
            if ((final_value == 1) and node_topoType == 'w'){
                png_ptr[i][4 * j] = f_R;
                png_ptr[i][4 * j + 1] = f_G;
                png_ptr[i][4 * j + 2] = f_B;
                png_ptr[i][4 * j + 3] = f_AD * final_value * 64 - 1;
            }
        }
    }
}

void dis_tileRender_classify0(int z, int x, int y, vector<string> attri, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, polygonNode *tile_node, polygonNode *RNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    vector<float> R0, R1, R2;
    for (int k = 0; k < width.size(); k++){
        R0.push_back(width[k] * Rz);
        R1.push_back(width[k] * Rz - sqrt(2) * Rz / 4.0);
        R2.push_back(width[k] * Rz + sqrt(2) * Rz / 4.0);
    }
    double Rz4 = 0.25 * Rz;
    float max_width = *max_element(width.begin(), width.end());
    int level = (max_width == 0) ? 0 : ceil(max_width + sqrt(2) / 4.0 - 0.5);
    // calcalate the BBox of tile <z,x,y>
    float tile_Rz = 256 * Rz;
    float tile_xMin = x * tile_Rz - L;
    float tile_xMax = (x + 1) * tile_Rz - L;
    float tile_yMin = L - (y + 1) * tile_Rz;
    float tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            // determine whether the node correspongding to the pixel exists
            polygonNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
            if (pixel_node){
            //     if (pixel_node->topo_type == 'i'){
            //         int ix = attriClassify<string>(pixel_node->attri_s, attri);
            //         if (ix != attri.size()){
            //             png_ptr[i][4 * j] = R[ix];     
            //             png_ptr[i][4 * j + 1] = G[ix]; 
            //             png_ptr[i][4 * j + 2] = B[ix]; 
            //             png_ptr[i][4 * j + 3] = AD[ix] * 255;
            //             continue;
            //         }
                // }
            //     else{
            //         png_ptr[i][4 * j] = 255;     
            //         png_ptr[i][4 * j + 1] = 255; 
            //         png_ptr[i][4 * j + 2] = 0; 
            //         png_ptr[i][4 * j + 3] = 255;
            //         continue;
            //     }
            }

            // short final_value = 0;
            // int final_ix;
            // short pix_value = 0;
            // int attri_ix;
            // for (int l = 1; l <= level; l++){
            //     for (int m = -l; m <= l; m++){
            //         for (int n = -l; n <= l; n++){
            //             if (abs(m) == l || abs(n) == l){
            //                 double nearPix_x = pix_x + m * Rz;
            //                 double nearPix_y = pix_y + n * Rz;
            //                 polygonNode *node_floor = locPixNode(z, tile_box, nearPix_x, nearPix_y, tile_node, RNode);
            //                 if (node_floor != nullptr){
            //                     attri_ix = attriClassify<string>(node_floor->attri_s, attri);
            //                     if (attri_ix == attri.size()){
            //                         continue;
            //                     }
                                
            //                     double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, nearPix_x, nearPix_y, Rz / 4);
            //                     double dist = nearInfo[0];
            //                     if (dist < R1[attri_ix]){
            //                         final_value = 4;
            //                         final_ix = attri_ix;
            //                         break;
            //                     }
            //                     else if (dist < R2[attri_ix]){  
            //                         double node_x = nearInfo[1];
            //                         double node_y = nearInfo[2];
            //                         if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
            //                             pix_value += 1;
            //                         }
            //                         if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
            //                             pix_value += 1;
            //                         }
            //                         if (sqrt(pow(pix_x + Rz4 - node_x, 2) + pow(pix_y - Rz4 - node_y, 2)) < R0[attri_ix]){
            //                             pix_value += 1;
            //                         }  
            //                         if (sqrt(pow(pix_x - Rz4 - node_x, 2) + pow(pix_y + Rz4 - node_y, 2)) < R0[attri_ix]){
            //                             pix_value += 1; 
            //                         }

            //                         if (pix_value == 4){
            //                             final_value = 4;
            //                             final_ix = attri_ix;
            //                             break;
            //                         }
            //                         else if (pix_value > final_value){
            //                             final_value = pix_value;
            //                             final_ix = attri_ix;
            //                         }
            //                         pix_value = 0;
            //                     }
            //                 }
            //             }
            //         }
            //         if (final_value == 4)
            //             break;
            //     }
            //     if (final_value == 4)
            //         break;
            // }
            // // generate pixel value
            // if (final_value > 0 and final_ix != attri.size()){
            //     png_ptr[i][4 * j] = R[final_ix];
            //     png_ptr[i][4 * j + 1] = G[final_ix];
            //     png_ptr[i][4 * j + 2] = B[final_ix];
            //     png_ptr[i][4 * j + 3] = AD[final_ix] * final_value * 64 - 1;
            // }
            // else{
            //     png_ptr[i][4 * j] = 0;
            //     png_ptr[i][4 * j + 1] = 0;
            //     png_ptr[i][4 * j + 2] = 0;
            //     png_ptr[i][4 * j + 3] = 0;
            // }
        }
    }
}

void dis_tileRender_classify(int z, int x, int y, vector<string> attri, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string field_name, polygonNode *tile_node, polygonNode *RNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    vector<float> R0, R1, R2;
    for (int k = 0; k < width.size(); k++){
        R0.push_back(width[k] * Rz);
        R1.push_back(width[k] * Rz - sqrt(2) * Rz / 4.0);
        R2.push_back(width[k] * Rz + sqrt(2) * Rz / 4.0);
    }
    double Rz4 = 0.25 * Rz;
    float max_width = *max_element(width.begin(), width.end());
    int level = (max_width == 0) ? 0 : ceil(max_width + sqrt(2) / 4.0 - 0.5);
    // calcalate the BBox of tile <z,x,y>
    float tile_Rz = 256 * Rz;
    float tile_xMin = x * tile_Rz - L;
    float tile_xMax = (x + 1) * tile_Rz - L;
    float tile_yMin = L - (y + 1) * tile_Rz;
    float tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            // determine whether the node correspongding to the pixel exists
            polygonNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
            if (pixel_node){
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
                        if (abs(m) == l || abs(n) == l){
                            double nearPix_x = pix_x + m * Rz;
                            double nearPix_y = pix_y + n * Rz;
                            polygonNode *node_floor = locPixNode(z, tile_box, nearPix_x, nearPix_y, tile_node, RNode);
                            if (node_floor != nullptr){
                                attri_ix = attriClassify(node_floor->attri_map[field_name], attri);
                                if (attri_ix == attri.size()){
                                    continue;
                                }
                                
                                double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, nearPix_x, nearPix_y, Rz / 4);
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

void dis_tileRender_range_classify(int z, int x, int y, double box[], vector<string> attri, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string field_name, polygonNode *tile_node, polygonNode *RNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    vector<float> R0, R1, R2;
    for (int k = 0; k < width.size(); k++){
        R0.push_back(width[k] * Rz);
        R1.push_back(width[k] * Rz - sqrt(2) * Rz / 4.0);
        R2.push_back(width[k] * Rz + sqrt(2) * Rz / 4.0);
    }
    double Rz4 = 0.25 * Rz;
    float max_width = *max_element(width.begin(), width.end());
    int level = (max_width == 0) ? 0 : ceil(max_width + sqrt(2) / 4.0 - 0.5);
    // calcalate the BBox of tile <z,x,y>
    float tile_Rz = 256 * Rz;
    float tile_xMin = x * tile_Rz - L;
    float tile_xMax = (x + 1) * tile_Rz - L;
    float tile_yMin = L - (y + 1) * tile_Rz;
    float tile_yMax = L - y * tile_Rz;
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
                polygonNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
                if (pixel_node){
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
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }

            short final_value = 0;
            int final_ix;
            short pix_value = 0;
            int attri_ix;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (abs(m) == l || abs(n) == l){
                            if (IfOverlapMBB(pix_xMin + m * Rz, pix_yMin + n * Rz, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                                double nearPix_x = pix_x + m * Rz;
                                double nearPix_y = pix_y + n * Rz;
                                polygonNode *node_floor = locPixNode(z, tile_box, nearPix_x, nearPix_y, tile_node, RNode);
                                if (node_floor != nullptr){
                                    attri_ix = attriClassify(node_floor->attri_map[field_name], attri);
                                    if (attri_ix == attri.size()){
                                        continue;
                                    }
                                    
                                    double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, nearPix_x, nearPix_y, Rz / 4);
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

void dis_tileRender_rule(int z, int x, int y, vector<string> filters, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string field_name, polygonNode *tile_node, polygonNode *RNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    vector<float> R0, R1, R2;
    for (int k = 0; k < width.size(); k++){
        R0.push_back(width[k] * Rz);
        R1.push_back(width[k] * Rz - sqrt(2) * Rz / 4.0);
        R2.push_back(width[k] * Rz + sqrt(2) * Rz / 4.0);
    }
    double Rz4 = 0.25 * Rz;
    float max_width = *max_element(width.begin(), width.end());
    int level = (max_width == 0) ? 0 : ceil(max_width + sqrt(2) / 4.0 - 0.5);
    // calcalate the BBox of tile <z,x,y>
    float tile_Rz = 256 * Rz;
    float tile_xMin = x * tile_Rz - L;
    float tile_xMax = (x + 1) * tile_Rz - L;
    float tile_yMin = L - (y + 1) * tile_Rz;
    float tile_yMax = L - y * tile_Rz;
    double tile_box[] = {tile_xMin, tile_yMin, tile_xMax, tile_yMax};
    // calcalate each pixel value
    # pragma omp parallel for num_threads(4)
    for (int i = 0; i < TILE_SIZE; i++){
        for (int j = 0; j < TILE_SIZE; j++){
            // calculate center coords of each pixel
            double pix_x = (256 * x + j + 0.5) * Rz - L;
            double pix_y = L - (256 * y + i + 0.5) * Rz;
            // determine whether the node correspongding to the pixel exists
            polygonNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
            if (pixel_node){
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
                        if (abs(m) == l || abs(n) == l){
                            double nearPix_x = pix_x + m * Rz;
                            double nearPix_y = pix_y + n * Rz;
                            polygonNode *node_floor = locPixNode(z, tile_box, nearPix_x, nearPix_y, tile_node, RNode);
                            if (node_floor != nullptr){
                                attri_ix = attriRule(node_floor->attri_map[field_name], filters);
                                if (attri_ix == filters.size()){
                                    continue;
                                }
                                
                                double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, nearPix_x, nearPix_y, Rz / 4);
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

void dis_tileRender_range_rule(int z, int x, int y, double box[], vector<string> filters, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string field_name, polygonNode *tile_node, polygonNode *RNode, png_bytep *png_ptr){
    double Rz = L / (128 << z);
    vector<float> R0, R1, R2;
    for (int k = 0; k < width.size(); k++){
        R0.push_back(width[k] * Rz);
        R1.push_back(width[k] * Rz - sqrt(2) * Rz / 4.0);
        R2.push_back(width[k] * Rz + sqrt(2) * Rz / 4.0);
    }
    double Rz4 = 0.25 * Rz;
    float max_width = *max_element(width.begin(), width.end());
    int level = (max_width == 0) ? 0 : ceil(max_width + sqrt(2) / 4.0 - 0.5);
    // calcalate the BBox of tile <z,x,y>
    float tile_Rz = 256 * Rz;
    float tile_xMin = x * tile_Rz - L;
    float tile_xMax = (x + 1) * tile_Rz - L;
    float tile_yMin = L - (y + 1) * tile_Rz;
    float tile_yMax = L - y * tile_Rz;
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
                polygonNode *pixel_node = nodeIfExist_attri(z, tile_box, pix_x, pix_y, tile_node);
                if (pixel_node){
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
            else{
                png_ptr[i][4 * j] = 0;
                png_ptr[i][4 * j + 1] = 0;
                png_ptr[i][4 * j + 2] = 0;
                png_ptr[i][4 * j + 3] = 0;
            }

            short final_value = 0;
            int final_ix;
            short pix_value = 0;
            int attri_ix;
            for (int l = 1; l <= level; l++){
                for (int m = -l; m <= l; m++){
                    for (int n = -l; n <= l; n++){
                        if (abs(m) == l || abs(n) == l){
                            if (IfOverlapMBB(pix_xMin + m * Rz, pix_yMin + n * Rz, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                                double nearPix_x = pix_x + m * Rz;
                                double nearPix_y = pix_y + n * Rz;
                                polygonNode *node_floor = locPixNode(z, tile_box, nearPix_x, nearPix_y, tile_node, RNode);
                                if (node_floor != nullptr){
                                    attri_ix = attriRule(node_floor->attri_map[field_name], filters);
                                    if (attri_ix == filters.size()){
                                        continue;
                                    }
                                    
                                    double *nearInfo = getNearestCNode(node_floor, pix_x, pix_y, nearPix_x, nearPix_y, Rz / 4);
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

// data-driven tile plotting
string data_tileRender_a_single(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, int R, int G, int B, float AD, float width, int f_R, int f_G, int f_B, float f_AD, string data_string, string edge_data_string){
    Map map(256, 256);
    map.set_srs(srs_merc);

    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt\n" + data_string;
    datasource_ptr polygon_ds = datasource_cache::instance().create(params);
    
    feature_type_style polygon_style;
    {
        rule r;
        {
            polygon_symbolizer poly_sym;
            put(poly_sym, keys::fill, color(f_R, f_G, f_B));
            put(poly_sym, keys::opacity, f_AD);
            put(poly_sym, keys::fill_opacity, f_AD);
            put(poly_sym, keys::gamma, 1.0);
            r.append(std::move(poly_sym));
        }
        polygon_style.add_rule(std::move(r));
    }
    map.insert_style("polygon_style", std::move(polygon_style));

    {
        layer lyr("polygon");
        lyr.set_datasource(polygon_ds);
        lyr.add_style("polygon_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    parameters edge_params;
    edge_params["type"] = "csv";
    edge_params["inline"] = "wkt\n" + edge_data_string;
    datasource_ptr line_ds = datasource_cache::instance().create(edge_params);

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

string data_tileRender_a_range_single(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, double box[], int R, int G, int B, float AD, float width, int f_R, int f_G, int f_B, float f_AD, string data_string, string edge_data_string){
    Map map(256, 256);
    map.set_srs(srs_merc);

    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt\n" + data_string;
    datasource_ptr polygon_ds = datasource_cache::instance().create(params);
    
    feature_type_style polygon_style;
    {
        rule r;
        {
            polygon_symbolizer poly_sym;
            put(poly_sym, keys::fill, color(f_R, f_G, f_B));
            put(poly_sym, keys::opacity, f_AD);
            put(poly_sym, keys::fill_opacity, f_AD);
            put(poly_sym, keys::gamma, 1.0);
            r.append(std::move(poly_sym));
        }
        polygon_style.add_rule(std::move(r));
    }
    map.insert_style("polygon_style", std::move(polygon_style));

    {
        layer lyr("polygon");
        lyr.set_datasource(polygon_ds);
        lyr.add_style("polygon_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    parameters edge_params;
    edge_params["type"] = "csv";
    edge_params["inline"] = "wkt\n" + edge_data_string;
    datasource_ptr line_ds = datasource_cache::instance().create(edge_params);

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

string data_tileRender_a_classify(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, string attri_name, string attri_type, vector<string> attri, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string data_string, string edge_data_string){
    Map map(256, 256);
    map.set_srs(srs_merc);

    // polygon style
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt," + attri_name + "\n" + data_string;
    datasource_ptr polygon_ds = datasource_cache::instance().create(params);
    
    feature_type_style polygon_style;
    polygon_style.reserve(attri.size());
    for (int i = 0; i < attri.size(); i++){
        rule r;
        string style_describer;

        if (attri_type == "string" or attri_type == "int"){
            style_describer = "[" + attri_name + "] = " + attri[i];
        }
        else if (attri_type == "real" or attri_type == "double"){
            style_describer = "[" + attri_name + "] = \'" + attri[i] + "\'";
        }

        r.set_filter(parse_expression(style_describer));
        {
            polygon_symbolizer poly_sym;
            put(poly_sym, keys::fill, color(R[i], G[i], B[i]));
            put(poly_sym, keys::opacity, AD[i]);
            put(poly_sym, keys::fill_opacity, AD[i]);
            put(poly_sym, keys::gamma, 1.0);
            r.append(std::move(poly_sym));
        }
        polygon_style.add_rule(std::move(r));
    }
    map.insert_style("polygon_style", std::move(polygon_style));

    {
        layer lyr("polygon");
        lyr.set_datasource(polygon_ds);
        lyr.add_style("polygon_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    // edge style
    // parameters edge_params;
    // edge_params["type"] = "csv";
    // edge_params["inline"] = "wkt\n" + edge_data_string;
    // datasource_ptr line_ds = datasource_cache::instance().create(edge_params);

    // feature_type_style line_style;
    // {
    //     rule r;
    //     {
    //         line_symbolizer line_sym;
    //         put(line_sym, keys::stroke, color(0, 0, 0));
    //         put(line_sym, keys::stroke_width, 1.0);
    //         put(line_sym, keys::stroke_linecap, ROUND_CAP);
    //         put(line_sym, keys::stroke_linejoin, ROUND_JOIN);
    //         put(line_sym, keys::stroke_opacity, 1.0);
    //         put(line_sym, keys::smooth, 1.0);
    //         r.append(std::move(line_sym));
    //     }
    //     line_style.add_rule(std::move(r));
    // }
    // map.insert_style("line_style", std::move(line_style));

    // {
    //     layer lyr("line");
    //     lyr.set_datasource(line_ds);
    //     lyr.add_style("line_style");
    //     lyr.set_srs(srs_merc);
    //     map.add_layer(lyr);
    // }

    map.zoom_to_box(box2d<double>(tile_xMin, tile_yMin, tile_xMax, tile_yMax));
    image_rgba8 im(256, 256);
    agg_renderer<image_rgba8> ren(map, im);
    ren.apply();
    ostringstream im_os;
    save_to_stream(im, im_os, "png");
    return im_os.str();
}

string data_tileRender_a_range_classify(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, double box[], string attri_name, string attri_type, vector<string> attri, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string data_string, string edge_data_string){
    Map map(256, 256);
    map.set_srs(srs_merc);

    // polygon style
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt," + attri_name + "\n" + data_string;
    datasource_ptr polygon_ds = datasource_cache::instance().create(params);
    
    feature_type_style polygon_style;
    polygon_style.reserve(attri.size());
    for (int i = 0; i < attri.size(); i++){
        rule r;
        string style_describer;

        if (attri_type == "string" or attri_type == "int"){
            style_describer = "[" + attri_name + "] = " + attri[i];
        }
        else if (attri_type == "real" or attri_type == "double"){
            style_describer = "[" + attri_name + "] = \'" + attri[i] + "\'";
        }

        r.set_filter(parse_expression(style_describer));
        {
            polygon_symbolizer poly_sym;
            put(poly_sym, keys::fill, color(R[i], G[i], B[i]));
            put(poly_sym, keys::opacity, AD[i]);
            put(poly_sym, keys::fill_opacity, AD[i]);
            put(poly_sym, keys::gamma, 1.0);
            r.append(std::move(poly_sym));
        }
        polygon_style.add_rule(std::move(r));
    }
    map.insert_style("polygon_style", std::move(polygon_style));

    {
        layer lyr("polygon");
        lyr.set_datasource(polygon_ds);
        lyr.add_style("polygon_style");
        lyr.set_srs(srs_merc);
        map.add_layer(lyr);
    }

    // edge style
    // parameters edge_params;
    // edge_params["type"] = "csv";
    // edge_params["inline"] = "wkt\n" + edge_data_string;
    // datasource_ptr line_ds = datasource_cache::instance().create(edge_params);

    // feature_type_style line_style;
    // {
    //     rule r;
    //     {
    //         line_symbolizer line_sym;
    //         put(line_sym, keys::stroke, color(0, 0, 0));
    //         put(line_sym, keys::stroke_width, 1.0);
    //         put(line_sym, keys::stroke_linecap, ROUND_CAP);
    //         put(line_sym, keys::stroke_linejoin, ROUND_JOIN);
    //         put(line_sym, keys::stroke_opacity, 1.0);
    //         put(line_sym, keys::smooth, 1.0);
    //         r.append(std::move(line_sym));
    //     }
    //     line_style.add_rule(std::move(r));
    // }
    // map.insert_style("line_style", std::move(line_style));

    // {
    //     layer lyr("line");
    //     lyr.set_datasource(line_ds);
    //     lyr.add_style("line_style");
    //     lyr.set_srs(srs_merc);
    //     map.add_layer(lyr);
    // }

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

string data_tileRender_a_rule(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, string attri_name, string attri_type, vector<string> filters, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string data_string){
    Map map(256, 256);
    map.set_srs(srs_merc);

    // polygon style
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt," + attri_name + "\n" + data_string;
    datasource_ptr polygon_ds = datasource_cache::instance().create(params);
    
    feature_type_style polygon_style;
    polygon_style.reserve(filters.size());
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
            polygon_symbolizer poly_sym;
            put(poly_sym, keys::fill, color(R[i], G[i], B[i]));
            put(poly_sym, keys::opacity, AD[i]);
            put(poly_sym, keys::fill_opacity, AD[i]);
            put(poly_sym, keys::gamma, 1.0);
            r.append(std::move(poly_sym));
        }
        polygon_style.add_rule(std::move(r));
    }
    map.insert_style("polygon_style", std::move(polygon_style));

    {
        layer lyr("polygon");
        lyr.set_datasource(polygon_ds);
        lyr.add_style("polygon_style");
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

string data_tileRender_a_range_rule(double tile_xMin, double tile_yMin, double tile_xMax, double tile_yMax, double box[], string attri_name, string attri_type, vector<string> filters, vector<int> R, vector<int> G, vector<int> B, vector<float> AD, vector<float> width, string data_string){
    Map map(256, 256);
    map.set_srs(srs_merc);

    // polygon style
    parameters params;
    params["type"] = "csv";
    params["inline"] = "wkt," + attri_name + "\n" + data_string;
    datasource_ptr polygon_ds = datasource_cache::instance().create(params);
    
    feature_type_style polygon_style;
    polygon_style.reserve(filters.size());
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
            polygon_symbolizer poly_sym;
            put(poly_sym, keys::fill, color(R[i], G[i], B[i]));
            put(poly_sym, keys::opacity, AD[i]);
            put(poly_sym, keys::fill_opacity, AD[i]);
            put(poly_sym, keys::gamma, 1.0);
            r.append(std::move(poly_sym));
        }
        polygon_style.add_rule(std::move(r));
    }
    map.insert_style("polygon_style", std::move(polygon_style));

    {
        layer lyr("polygon");
        lyr.set_datasource(polygon_ds);
        lyr.add_style("polygon_style");
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



// dis-driven visualization
string disDrivenPlot(int z, int x, int y, string style_describer, polygonNode *RNode){
    png_bytep * tile_ptr = (png_bytep *) malloc(256 * sizeof(png_bytep));
    for (int n = 0; n < TILE_SIZE; n++)
        tile_ptr[n] = (png_bytep) malloc(1024);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr); 
    FILE *temp_png = tmpfile();
    png_init_io(png_ptr, temp_png); 
    png_set_IHDR(png_ptr, info_ptr, TILE_SIZE, TILE_SIZE, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    polygonNode *tile_node = locTileNode(z, x, y, RNode);
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
            int f_R = stoi(style_list[5]);
            int f_G = stoi(style_list[6]);
            int f_B = stoi(style_list[7]);
            float f_AD = stof(style_list[8]);

            if (tile_node->topo_type == 'w'){
                # pragma omp parallel for num_threads(4)
                for (int i = 0; i < 256; i++){
                    for (int j = 0; j < 256; j++){
                        tile_ptr[i][4 * j] = f_R;
                        tile_ptr[i][4 * j + 1] = f_G;
                        tile_ptr[i][4 * j + 2] = f_B;
                        tile_ptr[i][4 * j + 3] = f_AD * 256 - 1;
                    }
                }
            }
            else{
                dis_tileRender_single(z, x, y, R, G, B, AD, width, f_R, f_G, f_B, f_AD, tile_node, RNode, tile_ptr);
            }
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
            
            if (tile_node->topo_type == 'w'){
                int ix = attriClassify(tile_node->attri_map[attri_name], attri);
                if (ix != attri.size()){
                    # pragma omp parallel for num_threads(4)
                    for (int i = 0; i < 256; i++){
                        for (int j = 0; j < 256; j++){
                            tile_ptr[i][4 * j] = R[ix];
                            tile_ptr[i][4 * j + 1] = G[ix];
                            tile_ptr[i][4 * j + 2] = B[ix];
                            tile_ptr[i][4 * j + 3] = AD[ix] * 256 - 1;
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
            }
            else{
                dis_tileRender_classify(z, x, y, attri, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
            }
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

            if (tile_node->topo_type == 'w'){
                int ix = attriRule(tile_node->attri_map[attri_name], filters);
                if (ix != filters.size()){
                    # pragma omp parallel for num_threads(4)
                    for (int i = 0; i < 256; i++){
                        for (int j = 0; j < 256; j++){
                            tile_ptr[i][4 * j] = R[ix];
                            tile_ptr[i][4 * j + 1] = G[ix];
                            tile_ptr[i][4 * j + 2] = B[ix];
                            tile_ptr[i][4 * j + 3] = AD[ix] * 256 - 1;
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
            }
            else{
                dis_tileRender_rule(z, x, y, filters, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
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

// dis-driven visualization
string disDrivenPlot_range(int z, int x, int y, double box[], string style_describer, polygonNode *RNode){
    png_bytep * tile_ptr = (png_bytep *) malloc(256 * sizeof(png_bytep));
    for (int n = 0; n < TILE_SIZE; n++)
        tile_ptr[n] = (png_bytep) malloc(1024);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr); 
    FILE *temp_png = tmpfile();
    png_init_io(png_ptr, temp_png); 
    png_set_IHDR(png_ptr, info_ptr, TILE_SIZE, TILE_SIZE, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    polygonNode *tile_node = locTileNode(z, x, y, RNode);
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
            int f_R = stoi(style_list[5]);
            int f_G = stoi(style_list[6]);
            int f_B = stoi(style_list[7]);
            float f_AD = stof(style_list[8]);

            if (tile_node->topo_type == 'w'){
                double Rz = L / (128 << z);
                # pragma omp parallel for num_threads(4)
                for (int i = 0; i < 256; i++){
                    for (int j = 0; j < 256; j++){
                        double pix_xMin = (256 * x + j) * Rz - L;
                        double pix_yMin = L - (256 * y + i + 1) * Rz;
                        if (IfOverlapMBB(pix_xMin, pix_yMin, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                            tile_ptr[i][4 * j] = f_R;
                            tile_ptr[i][4 * j + 1] = f_G;
                            tile_ptr[i][4 * j + 2] = f_B;
                            tile_ptr[i][4 * j + 3] = f_AD * 256 - 1;
                        }
                        else{
                            tile_ptr[i][4 * j] = 0;
                            tile_ptr[i][4 * j + 1] = 0;
                            tile_ptr[i][4 * j + 2] = 0;
                            tile_ptr[i][4 * j + 3] = 0;
                        }
                    }
                }
            }
            else{
                if (IfContain(box_xMin, box_yMin, box_xMax, box_yMax, tile_xMin, tile_yMin, tile_xMin+tile_Rz, tile_yMin+tile_Rz)){
                    dis_tileRender_single(z, x, y, R, G, B, AD, width, f_R, f_G, f_B, f_AD, tile_node, RNode, tile_ptr);
                }
                else{
                    dis_tileRender_range_single(z, x, y, box, R, G, B, AD, width, f_R, f_G, f_B, f_AD, tile_node, RNode, tile_ptr);
                }
            }
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
            
            if (tile_node->topo_type == 'w'){
                int ix = attriClassify(tile_node->attri_map[attri_name], attri);
                if (ix != attri.size()){
                    double Rz = L / (128 << z);
                    # pragma omp parallel for num_threads(4)
                    for (int i = 0; i < 256; i++){
                        for (int j = 0; j < 256; j++){
                            double pix_xMin = (256 * x + j) * Rz - L;
                            double pix_yMin = L - (256 * y + i + 1) * Rz;
                            if (IfOverlapMBB(pix_xMin, pix_yMin, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                                tile_ptr[i][4 * j] = R[ix];
                                tile_ptr[i][4 * j + 1] = G[ix];
                                tile_ptr[i][4 * j + 2] = B[ix];
                                tile_ptr[i][4 * j + 3] = AD[ix] * 256 - 1;
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
            }
            else{
                if (IfContain(box_xMin, box_yMin, box_xMax, box_yMax, tile_xMin, tile_yMin, tile_xMin+tile_Rz, tile_yMin+tile_Rz)){
                    dis_tileRender_classify(z, x, y, attri, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
                }
                else{
                    dis_tileRender_range_classify(z, x, y, box, attri, R, G, B, AD, width, attri_name, tile_node, RNode, tile_ptr);
                }
            }
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

            if (tile_node->topo_type == 'w'){
                int ix = attriRule(tile_node->attri_map[attri_name], filters);
                if (ix != filters.size()){
                    double Rz = L / (128 << z);
                    # pragma omp parallel for num_threads(4)
                    for (int i = 0; i < 256; i++){
                        for (int j = 0; j < 256; j++){
                            double pix_xMin = (256 * x + j) * Rz - L;
                            double pix_yMin = L - (256 * y + i + 1) * Rz;
                            if (IfOverlapMBB(pix_xMin, pix_yMin, Rz, box_xMin, box_yMin, box_xMax, box_yMax)){
                                tile_ptr[i][4 * j] = R[ix];
                                tile_ptr[i][4 * j + 1] = G[ix];
                                tile_ptr[i][4 * j + 2] = B[ix];
                                tile_ptr[i][4 * j + 3] = AD[ix] * 256 - 1;
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
            }
            else{
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
string dataDrivenPlot(int z, int x, int y, string style_describer, polygonNode *RNode, short switch_level){
    string img_str;
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    float tile_xCenter = tile_xMin + (tile_xMax - tile_xMin) / 2;
    float tile_yCenter = tile_yMin + (tile_yMax - tile_yMin) / 2;
    // determine whether the tile needs to be drawn
    pair<polygonNode*, polygonNode*> node_ptrs = locRtreeNode(z, tile_xCenter, tile_yCenter, RNode, switch_level);
    polygonNode *tile_node = node_ptrs.first;
    // polygonNode *rtree_node = locRtreeNode(z, tile_xCenter, tile_yCenter, RNode, switch_level);
    // non-blank tile
    if (tile_node != nullptr){
        string data_string = "";
        string edge_data_string = "";
        vector<string> attri_list;
        boost::split(attri_list, style_describer, boost::is_any_of( "-" ), boost::token_compress_on);
        string style_type = attri_list[0];

        // get the features by spatial range retrieve
        vector<POLYGON_ITEM> polys;
        Box box(POINT(tile_xMin, tile_yMin), POINT(tile_xMax, tile_yMax));
        if (tile_node->topo_type == 'i'){
            polygonNode *rtree_node = node_ptrs.second;
            rtree_node->Rtree.query(bgi::intersects(box), back_inserter(polys));


            // single visualization
            if (style_type == "single"){
                vector<string> style_list;
                boost::split(style_list, attri_list[1], boost::is_any_of( "," ), boost::token_compress_on);
                int R = stoi(style_list[0]);
                int G = stoi(style_list[1]);
                int B = stoi(style_list[2]);
                float AD = stof(style_list[3]);
                float width = stof(style_list[4]);
                int f_R = stoi(style_list[5]);
                int f_G = stoi(style_list[6]);
                int f_B = stoi(style_list[7]);
                float f_AD = stof(style_list[8]);
                
                // sort(begin(polys), end(polys), [](const POLYGON& p1, const POLYGON& p2) {return bg::area(get<0>(p1)) > bg::area(get<0>(p2));});
                for (int i = 0; i < polys.size(); i++){
                    RING ring_cord = get<1>(polys[i]);
                    data_string += "\"POLYGON((";
                    for (int j = 0; j < ring_cord.size()-2; j+=2){
                        data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                        if (j <= ring_cord.size()-4){
                            edge_data_string += "\"LINESTRING(";
                            edge_data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                            edge_data_string += (to_string(ring_cord[j+2]) + " " + to_string(ring_cord[j+3]) + ")\"\n");
                        }
                    }
                    data_string += (to_string(ring_cord[0]) + " " + to_string(ring_cord[1]) + "))\"\n");
                }
                // plot features
                img_str = data_tileRender_a_single(tile_xMin, tile_yMin, tile_xMax, tile_yMax, R, G, B, AD, width, f_R, f_G, f_B, f_AD, data_string, edge_data_string);
            }
            else{
                vector<int> R, G, B;
                vector<float> AD, width;
                string attri_name = attri_list[1];
                string attri_type = attri_list[2];

                // sort(begin(polys), end(polys), [](const POLYGON& p1, const POLYGON& p2) {return bg::area(get<0>(p1)) > bg::area(get<0>(p2));});
                for (int i = 0; i < polys.size(); i++){
                    RING ring_cord = get<1>(polys[i]);
                    polygonAttri attri = get<2>(polys[i]);
                    data_string += "\"POLYGON((";
                    for (int j = 0; j < ring_cord.size()-2; j+=2){
                        data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                        if (j <= ring_cord.size()-4){
                            edge_data_string += "\"LINESTRING(";
                            edge_data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                            edge_data_string += (to_string(ring_cord[j+2]) + " " + to_string(ring_cord[j+3]) + ")\"\n");
                        }
                    }
                    data_string += (to_string(ring_cord[0]) + " " + to_string(ring_cord[1]) + "))\",");
                    // if (attri_type == "string"){
                    //     data_string += (attri.field_s + "\n");
                    // }
                    // else if (attri_type == "int"){
                    //     data_string += (to_string(attri.field_i) + "\n");
                    // }
                    // else if (attri_type == "double"){
                    //     data_string += (to_string(attri.field_d) + "\n");
                    // }
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
                    img_str = data_tileRender_a_classify(tile_xMin, tile_yMin, tile_xMax, tile_yMax, attri_name, attri_type, attri, R, G, B, AD, width, data_string, edge_data_string);
                }
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
                    img_str = data_tileRender_a_rule(tile_xMin, tile_yMin, tile_xMax, tile_yMax, attri_name, attri_type, rules, R, G, B, AD, width, data_string);
                }
            }
        }
        else if (tile_node->topo_type == 'w'){
            // tile_node->Rtree.query(bgi::intersects(box), back_inserter(polys));
            int f_R = 0;
            int f_G = 0;
            int f_B = 0;
            float f_AD = 0;
            if (style_type == "single"){
                vector<string> style_list;
                boost::split(style_list, attri_list[1], boost::is_any_of( "," ), boost::token_compress_on);
                f_R = stoi(style_list[5]);
                f_G = stoi(style_list[6]);
                f_B = stoi(style_list[7]);
                f_AD = stof(style_list[8]);
            }
            else if (style_type == "classified"){
                string attri_name = attri_list[1];
                string field_item = tile_node->attri_map[attri_name].begin()->first;

                for (int i = 3; i < attri_list.size(); i++){
                    vector<string> style_list;
                    boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                    if (field_item == style_list[0]){
                        f_R = stoi(style_list[1]);
                        f_G = stoi(style_list[2]);
                        f_B = stoi(style_list[3]);
                        f_AD = stof(style_list[4]);
                        break;
                    }
                }
            }
            else if (style_type == "ruled"){
                vector <string> rules, R, G, B, AD;
                for (int i = 3; i < attri_list.size(); i++){
                    vector<string> style_list;
                    boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                    rules.push_back(style_list[0]);
                    R.push_back(style_list[1]);
                    G.push_back(style_list[2]);
                    B.push_back(style_list[3]);
                    AD.push_back(style_list[4]);
                }
                
                string attri_name = attri_list[1];

                int ix = attriRule(tile_node->attri_map[attri_name], rules);
                if (ix != rules.size()){
                    f_R = stoi(R[ix]);
                    f_G = stoi(G[ix]);
                    f_B = stoi(B[ix]);
                    f_AD = stof(AD[ix]);
                }
            }

            Map map(256, 256);
            const color back_color(f_R, f_G, f_B, 255 * f_AD);
            map.set_background(back_color);
            image_rgba8 im(256, 256);
            agg_renderer<image_rgba8> ren(map, im);
            ren.apply();
            ostringstream im_os;
            save_to_stream(im, im_os, "png");
            img_str = im_os.str();
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

string dataDrivenPlot_range(int z, int x, int y, double range_box[], string style_describer, polygonNode *RNode, short switch_level){
    string img_str;
    double tile_Rz = L / pow(2, z-1);
    double tile_xMin = x * tile_Rz - L;
    double tile_xMax = (x + 1) * tile_Rz - L;
    double tile_yMin = L - (y + 1) * tile_Rz;
    double tile_yMax = L - y * tile_Rz;
    float tile_xCenter = tile_xMin + (tile_xMax - tile_xMin) / 2;
    float tile_yCenter = tile_yMin + (tile_yMax - tile_yMin) / 2;
    // determine whether the tile needs to be drawn
    pair<polygonNode*, polygonNode*> node_ptrs = locRtreeNode(z, tile_xCenter, tile_yCenter, RNode, switch_level);
    polygonNode *tile_node = node_ptrs.first;
    // polygonNode *rtree_node = locRtreeNode(z, tile_xCenter, tile_yCenter, RNode, switch_level);
    // non-blank tile
    if (tile_node != nullptr){
        double box_xMin = range_box[0];
        double box_yMin = range_box[1];
        double box_xMax = range_box[2];
        double box_yMax = range_box[3];
        
        string data_string = "";
        string edge_data_string = "";
        vector<string> attri_list;
        boost::split(attri_list, style_describer, boost::is_any_of( "-" ), boost::token_compress_on);
        string style_type = attri_list[0];

        // get the features by spatial range retrieve
        vector<POLYGON_ITEM> polys;
        Box box(POINT(tile_xMin, tile_yMin), POINT(tile_xMax, tile_yMax));
        if (tile_node->topo_type == 'i'){
            polygonNode *rtree_node = node_ptrs.second;
            rtree_node->Rtree.query(bgi::intersects(box), back_inserter(polys));


            // single visualization
            if (style_type == "single"){
                vector<string> style_list;
                boost::split(style_list, attri_list[1], boost::is_any_of( "," ), boost::token_compress_on);
                int R = stoi(style_list[0]);
                int G = stoi(style_list[1]);
                int B = stoi(style_list[2]);
                float AD = stof(style_list[3]);
                float width = stof(style_list[4]);
                int f_R = stoi(style_list[5]);
                int f_G = stoi(style_list[6]);
                int f_B = stoi(style_list[7]);
                float f_AD = stof(style_list[8]);
                
                // sort(begin(polys), end(polys), [](const POLYGON& p1, const POLYGON& p2) {return bg::area(get<0>(p1)) > bg::area(get<0>(p2));});
                for (int i = 0; i < polys.size(); i++){
                    RING ring_cord = get<1>(polys[i]);
                    data_string += "\"POLYGON((";
                    for (int j = 0; j < ring_cord.size()-2; j+=2){
                        data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                        if (j <= ring_cord.size()-4){
                            edge_data_string += "\"LINESTRING(";
                            edge_data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                            edge_data_string += (to_string(ring_cord[j+2]) + " " + to_string(ring_cord[j+3]) + ")\"\n");
                        }
                    }
                    data_string += (to_string(ring_cord[0]) + " " + to_string(ring_cord[1]) + "))\"\n");
                }

                // plot features
                if (IfContain(box_xMin, box_yMin, box_xMax, box_yMax, tile_xMin, tile_yMin, tile_xMax, tile_yMax)){
                    img_str = data_tileRender_a_single(tile_xMin, tile_yMin, tile_xMax, tile_yMax, R, G, B, AD, width, f_R, f_G, f_B, f_AD, data_string, edge_data_string);
                }
                else{
                    img_str = data_tileRender_a_range_single(tile_xMin, tile_yMin, tile_xMax, tile_yMax, range_box, R, G, B, AD, width, f_R, f_G, f_B, f_AD, data_string, edge_data_string);
                }
            }
            else{
                vector<int> R, G, B;
                vector<float> AD, width;
                string attri_name = attri_list[1];
                string attri_type = attri_list[2];

                // sort(begin(polys), end(polys), [](const POLYGON& p1, const POLYGON& p2) {return bg::area(get<0>(p1)) > bg::area(get<0>(p2));});
                for (int i = 0; i < polys.size(); i++){
                    RING ring_cord = get<1>(polys[i]);
                    polygonAttri attri = get<2>(polys[i]);
                    data_string += "\"POLYGON((";
                    for (int j = 0; j < ring_cord.size()-2; j+=2){
                        data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                        if (j <= ring_cord.size()-4){
                            edge_data_string += "\"LINESTRING(";
                            edge_data_string += (to_string(ring_cord[j]) + " " + to_string(ring_cord[j+1]) + ",");
                            edge_data_string += (to_string(ring_cord[j+2]) + " " + to_string(ring_cord[j+3]) + ")\"\n");
                        }
                    }
                    data_string += (to_string(ring_cord[0]) + " " + to_string(ring_cord[1]) + "))\",");
                    // if (attri_type == "string"){
                    //     data_string += (attri.field_s + "\n");
                    // }
                    // else if (attri_type == "int"){
                    //     data_string += (to_string(attri.field_i) + "\n");
                    // }
                    // else if (attri_type == "double"){
                    //     data_string += (to_string(attri.field_d) + "\n");
                    // }
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
                        img_str = data_tileRender_a_classify(tile_xMin, tile_yMin, tile_xMax, tile_yMax, attri_name, attri_type, attri, R, G, B, AD, width, data_string, edge_data_string);
                    }
                    else{
                        img_str = data_tileRender_a_range_classify(tile_xMin, tile_yMin, tile_xMax, tile_yMax, range_box, attri_name, attri_type, attri, R, G, B, AD, width, data_string, edge_data_string);
                    }
                }
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
                        img_str = data_tileRender_a_rule(tile_xMin, tile_yMin, tile_xMax, tile_yMax, attri_name, attri_type, rules, R, G, B, AD, width, data_string);
                    }
                    else{
                        img_str = data_tileRender_a_range_rule(tile_xMin, tile_yMin, tile_xMax, tile_yMax, range_box, attri_name, attri_type, rules, R, G, B, AD, width, data_string);
                    }
                }
            }
        }
        else if (tile_node->topo_type == 'w'){
            // tile_node->Rtree.query(bgi::intersects(box), back_inserter(polys));
            int f_R = 0;
            int f_G = 0;
            int f_B = 0;
            float f_AD = 0;
            if (style_type == "single"){
                vector<string> style_list;
                boost::split(style_list, attri_list[1], boost::is_any_of( "," ), boost::token_compress_on);
                f_R = stoi(style_list[5]);
                f_G = stoi(style_list[6]);
                f_B = stoi(style_list[7]);
                f_AD = stof(style_list[8]);
            }
            else if (style_type == "classified"){
                string attri_name = attri_list[1];
                string field_item = tile_node->attri_map[attri_name].begin()->first;

                for (int i = 3; i < attri_list.size(); i++){
                    vector<string> style_list;
                    boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                    if (field_item == style_list[0]){
                        f_R = stoi(style_list[1]);
                        f_G = stoi(style_list[2]);
                        f_B = stoi(style_list[3]);
                        f_AD = stof(style_list[4]);
                        break;
                    }
                }
            }
            else if (style_type == "ruled"){
                vector <string> rules, R, G, B, AD;
                for (int i = 3; i < attri_list.size(); i++){
                    vector<string> style_list;
                    boost::split(style_list, attri_list[i], boost::is_any_of( "," ), boost::token_compress_on);
                    rules.push_back(style_list[0]);
                    R.push_back(style_list[1]);
                    G.push_back(style_list[2]);
                    B.push_back(style_list[3]);
                    AD.push_back(style_list[4]);
                }
                
                string attri_name = attri_list[1];

                int ix = attriRule(tile_node->attri_map[attri_name], rules);
                if (ix != rules.size()){
                    f_R = stoi(R[ix]);
                    f_G = stoi(G[ix]);
                    f_B = stoi(B[ix]);
                    f_AD = stof(AD[ix]);
                }
            }

            Map map(256, 256);
            const color back_color(f_R, f_G, f_B, 255 * f_AD);
            map.set_background(back_color);
            image_rgba8 im(256, 256);
            agg_renderer<image_rgba8> ren(map, im);
            ren.apply();

            double pix_Rz = tile_Rz / 256;
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
            img_str = im_os.str();
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

string fieldStatistic(double box_xMin, double box_yMin, double box_xMax, double box_yMax, string style_describer, short switch_level, polygonNode *RNode){
    string field_statistic = "";
    map<string, double> count_list;
    vector<string> attri_list;
    vector<string> style_list;
    boost::split(attri_list, style_describer, boost::is_any_of( "-" ), boost::token_compress_on);
    string style_type = attri_list[0];

    if (style_type == "single"){
        polygonNode *min_WNode = RNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);

        if (min_WNode){
            min_WNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_list);
            if (count_list.size()){
                field_statistic = "Count:" + to_string(count_list["area"]);
            }
            else{
                field_statistic = "Count:0";
            }
        }
        else{
            field_statistic = "Count:0";
        }
    }
    else if (style_type == "classified"){
        string attri_name = attri_list[1];
        polygonNode *min_WNode = RNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
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
        polygonNode *min_WNode = RNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        if (min_WNode){
            vector<string> rules;
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



