#include <iostream>
#include "treeNode.hpp"
using namespace std;

#ifndef _POLYGONPLOT_HPP_
#define _POLYGONPLOT_HPP_


string disDrivenPlot(int z , int x, int y, string style_describer, polygonNode *RNode);
string disDrivenPlot_range(int z, int x, int y, double box[], string style_describer, polygonNode *RNode);

string dataDrivenPlot(int z , int x, int y, string style_describer, polygonNode *RNode, short switch_level);
string dataDrivenPlot_range(int z, int x, int y, double range_box[], string style_describer, polygonNode *RNode, short switch_level);

polygonNode * locTileNode(int z, int x, int y, polygonNode *seek);

string fieldStatistic(double box_xMin, double box_yMin, double box_xMax, double box_yMax, string style_describer, short switch_level, polygonNode *RNode);

#endif