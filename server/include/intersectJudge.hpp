#include <iostream>

#ifndef INTERSECTJUDGE_HPP_
#define INTERSECTJUDGE_HPP_

using namespace std;


bool IfContainPoint(double node_x, double node_y, double node_w, double obj_x, double obj_y);
bool IfContainMBB(double node_x, double node_y, double node_w, double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax);
bool IfContain(double c_xMin, double c_yMin, double c_xMax, double c_yMax, double w_xMin, double w_yMin, double w_xMax, double w_yMax);
bool IfOverlapMBB(double node_x, double node_y, double width, double obj_x_min, double obj_y_min, double obj_x_max, double obj_y_max);
bool IfInterctSeg(double node_x, double node_y, double node_w, double seg_xMin, double seg_yMin, double seg_xMax, double seg_yMax, string seg_type);
bool IfwithinBox(float x, float y, float box_xMin, float box_yMin, float box_xMax, float box_yMax);
bool IfwithinBox(double node_x, double node_y, double node_w, double box_xMin, double box_yMin, double box_xMax, double box_yMax);
bool IfwithinPolygon(float x, float y, float ring_x[], float ring_y[], int len);

#endif