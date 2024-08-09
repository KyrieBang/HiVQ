#include <iostream>
#include "treeNode.hpp"
using namespace std;

#ifndef _POINTPLOT_HPP_
#define _POINTPLOT_HPP_

typedef array<array<int, 3>, 256> RGB;


string disDrivenPlot(int z , int x, int y, string style_describer, pointNode *pointRNode);

string disDrivenPlot_range(int z, int x, int y, double box[], string style_describer, pointNode *pointRNode);

string dataDrivenPlot(int z , int x, int y, string style_describer, pointNode *pointRNode, short switch_level);

string dataDrivenPlot_range(int z, int x, int y, double range_box[], string style_describer, pointNode *pointRNode, short switch_level);

string fieldStatistic(double box_xMin, double box_yMin, double box_xMax, double box_yMax, string style_describer, short switch_level, pointNode *RNode);

pointNode * locTileNode(int z, int x, int y, pointNode *seek);

bool nodeIfExist(int z, double tile_box[], double pix_x, double pix_y, pointNode *seek);

bool locPixNode(int z, double tile_box[], double pix_x, double pix_y, pointNode *seek, pointNode *rootNode);

RGB initRGBTable(string gradient_type);

#endif