#include <iostream>
#include "treeNode.hpp"
using namespace std;

#ifndef _LINESTRINGPLOT_HPP_
#define _LINESTRINGPLOT_HPP_


string disDrivenPlot(int z , int x, int y, string style_describer, linestringNode *RNode);

string disDrivenPlot_range(int z, int x, int y, double box[], string style_describer, linestringNode *RNode);

string dataDrivenPlot(int z , int x, int y, string style_describer, linestringNode *RNode, short switch_level);

string dataDrivenPlot_range(int z, int x, int y, double range_box[], string style_describer, linestringNode *RNode, short switch_level);

bool nodeIfExist(int z, double tile_box[], double pix_x, double pix_y, linestringNode *seek);

linestringNode * locPixNode(int z, double tile_box[], double pix_x, double pix_y, linestringNode *seek, linestringNode *RNode);

double *getNearestCNode(linestringNode *seek, double pix_x, double pix_y, double near_x, double near_y, double Rz);

string fieldStatistic(double box_xMin, double box_yMin, double box_xMax, double box_yMax, string style_describer, short switch_level, linestringNode *RNode);

linestringNode * locTileNode(int z, int x, int y, linestringNode *seek);

int attriClassify(map<string, double> node_attri, vector<string> attri_list);

int attriRule(map<string, double> node_attri, vector<string> filter_list);

bool cmp_value(const pair<string, double> left, const pair<string, double> right);

#endif