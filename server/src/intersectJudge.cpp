#include "intersectJudge.hpp"
#include <iostream>
using namespace std;



// determine whether the BBox of the index node contains the obj.point
bool IfContainPoint(double node_x, double node_y, double node_w, double obj_x, double obj_y){
    if (obj_x >= node_x && 
        obj_x <= node_x + node_w && 
        obj_y >= node_y && 
        obj_y <= node_y + node_w){
        return true;
    }
    return false;
}

// determine whether the BBox of the index node contains the obj.MBB
bool IfContainMBB(double node_x, double node_y, double node_w, double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax){
    if (obj_xMin >= node_x && 
        obj_xMax <= node_x + node_w && 
        obj_yMin >= node_y && 
        obj_yMax <= node_y + node_w){
        return true;
    }
    return false;
}

bool IfContain(double c_xMin, double c_yMin, double c_xMax, double c_yMax, double w_xMin, double w_yMin, double w_xMax, double w_yMax){
    if (w_xMin >= c_xMin && w_xMax <= c_xMax && w_yMin >= c_yMin && w_yMax <= c_yMax){
        return true;
    }
    return false;
}

// determine whether the BBox of the index node overlaps the obj.MBB
bool IfOverlapMBB(double node_x, double node_y, double node_w, double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax){
    if (obj_xMin < node_x + node_w && 
        obj_xMax > node_x && 
        obj_yMin < node_y + node_w && 
        obj_yMax > node_y){
        return true;
    }
    return false;
}


// determine whether the BBox of the index node overlaps the obj.seg
bool IfInterctSeg(double node_x, double node_y, double node_w, double seg_xMin, double seg_yMin, double seg_xMax, double seg_yMax, string seg_type){
    if (seg_xMin < node_x + node_w && 
        seg_xMax > node_x && 
        seg_yMin < node_y + node_w && 
        seg_yMax > node_y){

        double center_x = node_x + node_w / 2;
        double center_y = node_y + node_w / 2;
        double seg_width = seg_xMax - seg_xMin;
        double seg_height = seg_yMax - seg_yMin;
        if (seg_type == "0"){
            double a = seg_width * (center_y - seg_yMin) - (center_x - seg_xMin) * seg_height;
            if (a == 0){
                return true;
            }
            else if (a > 0){
                double b = seg_width * (node_y - seg_yMin) - (node_x + node_w - seg_xMin) * seg_height;
                if (a * b < 0){
                    return true;
                }
            }
            else if (a < 0){
                double b = seg_width * (node_y + node_w - seg_yMin) - (node_x - seg_xMin) * seg_height;
                if (a * b < 0){
                    return true;
                }
            }
        }
        else if (seg_type == "1"){
            double a = seg_width * (center_y - seg_yMin - seg_height) + (center_x - seg_xMin) * seg_height;
            if (a == 0){
                return true;
            }
            else if (a > 0){
                double b = seg_width * (node_y - seg_yMin - seg_height) + (node_x - seg_xMin) * seg_height;
                if (a * b < 0){
                    return true;
                }
            }
            else if (a < 0){
                double b = seg_width * (node_y + node_w - seg_yMin - seg_height) + (node_x + node_w - seg_xMin) * seg_height;
                if (a * b < 0){
                    return true;
                }
            }

        }
        return false;
    }
    return false;
}

// determine whether the point within the box
bool IfwithinBox(float x, float y, float box_xMin, float box_yMin, float box_xMax, float box_yMax){
    if (x < box_xMin || x > box_xMax || y < box_yMin || y > box_yMax){
        return false;
    }
    return true;
}

bool IfwithinBox(double node_x, double node_y, double node_w, double box_xMin, double box_yMin, double box_xMax, double box_yMax){
    if (box_xMin <= node_x && 
        box_xMax >= node_x + node_w && 
        box_yMin <= node_y && 
        box_yMax >= node_y + node_w){
        return true;
    }
    return false;
}

// determine whether the point whtin the polygon
bool IfwithinPolygon(float x, float y, float ring_x[], float ring_y[], int len){
    int lcount = 0;
    float diff_x0 = ring_x[0] - x;
    float diff_y0 = ring_y[0] - y;
    for (int i = 1; i < len; i++){
        float diff_x1 = ring_x[i] - x;
        float diff_y1 = ring_y[i] - y;
        if (((diff_y0 <= 0) && (diff_y1 > 0)) || ((diff_y0 > 0) && (diff_y1 <= 0))){
            const float intersection = (diff_x1 * diff_y0 - diff_x0 * diff_y1) / (diff_y0 - diff_y1);
            if (intersection > 0){
                lcount++;
            }
        }
        diff_x0 = diff_x1;
        diff_y0 = diff_y1;
    }
   
    if (lcount % 2 == 1){
        return true;
    } 
    
    return false;
}