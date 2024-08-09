#include "treeNode.hpp"
#include "intersectJudge.hpp"
#include <iostream>
#include <fstream>
#include <sys/mman.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry/geometries/box.hpp>

#define L 20037508.3427892


void updateCount_ruled(map<string, int> fields, vector<string> rules, map<string, int> &count_map){
    int ix;
    for (auto rule_i: rules){
        for (auto field_i: fields){
            string field_i_key = field_i.first;
            int field_i_value = field_i.second;
            if ((ix = rule_i.find("!=")) != string::npos){
                if (field_i_key != rule_i.substr(ix+2)){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find("%3E=")) != string::npos){
                if (stod(field_i_key) >= stod(rule_i.substr(ix+4))){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find("%3C=")) != string::npos){
                if (stod(field_i_key) <= stod(rule_i.substr(ix+4))){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find('=')) != string::npos){
                if (field_i_key == rule_i.substr(ix+1)){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find("%3E")) != string::npos){
                if (stod(field_i_key) > stod(rule_i.substr(ix+3))){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find("%3C")) != string::npos){
                if (stod(field_i_key) < stod(rule_i.substr(ix+3))){
                    count_map[rule_i] += field_i_value;
                }
            }
        }
    }
}


void updateCount_ruled(map<string, double> fields, vector<string> rules, map<string, double> &count_map){
    int ix;
    for (auto rule_i: rules){
        for (auto field_i: fields){
            string field_i_key = field_i.first;
            double field_i_value = field_i.second;
            if ((ix = rule_i.find("!=")) != string::npos){
                if (field_i_key != rule_i.substr(ix+2)){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find("%3E=")) != string::npos){
                if (stod(field_i_key) >= stod(rule_i.substr(ix+4))){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find("%3C=")) != string::npos){
                if (stod(field_i_key) <= stod(rule_i.substr(ix+4))){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find('=')) != string::npos){
                if (field_i_key == rule_i.substr(ix+1)){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find("%3E")) != string::npos){
                if (stod(field_i_key) > stod(rule_i.substr(ix+3))){
                    count_map[rule_i] += field_i_value;
                }
            }
            else if ((ix = rule_i.find("%3C")) != string::npos){
                if (stod(field_i_key) < stod(rule_i.substr(ix+3))){
                    count_map[rule_i] += field_i_value;
                }
            }
        }
    }
}


void updateCount_withinNode_ruled(map<string, double> fields, double i_area, vector<string> rules, map<string, double> &count_map){
    int ix;
    for (auto rule_i: rules){
        for (auto field_i: fields){
            string field_i_key = field_i.first;
            if ((ix = rule_i.find("!=")) != string::npos){
                if (field_i_key != rule_i.substr(ix+2)){
                    count_map[rule_i] += i_area;
                }
            }
            else if ((ix = rule_i.find("%3E=")) != string::npos){
                if (stod(field_i_key) >= stod(rule_i.substr(ix+4))){
                    count_map[rule_i] += i_area;
                }
            }
            else if ((ix = rule_i.find("%3C=")) != string::npos){
                if (stod(field_i_key) <= stod(rule_i.substr(ix+4))){
                    count_map[rule_i] += i_area;
                }
            }
            else if ((ix = rule_i.find('=')) != string::npos){
                if (field_i_key == rule_i.substr(ix+1)){
                    count_map[rule_i] += i_area;
                }
            }
            else if ((ix = rule_i.find("%3E")) != string::npos){
                if (stod(field_i_key) > stod(rule_i.substr(ix+3))){
                    count_map[rule_i] += i_area;
                }
            }
            else if ((ix = rule_i.find("%3C")) != string::npos){
                if (stod(field_i_key) < stod(rule_i.substr(ix+3))){
                    count_map[rule_i] += i_area;
                }
            }
        }
    }
}



// ---------------------------------------------------------------------------------------------------------------------------------------------- //
// ***************************************************    point node functions    *************************************************************** //
// ---------------------------------------------------------------------------------------------------------------------------------------------- //
pointNode::pointNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, pointNode *_parent){
    this->node_type = _node_type;
    this->node_x = _node_x;
    this->node_y = _node_y;
    this->node_w = _node_w;
    this->node_level = _node_level;
    this->node_fcount = 0;
    this->parent = _parent;
    this->LUNode = nullptr;
    this->RUNode = nullptr;
    this->LBNode = nullptr;
    this->RBNode = nullptr;
}


pointNode::~pointNode()
{
    delete LUNode;
    delete RUNode;
    delete LBNode;
    delete RBNode;
    delete parent;
    LUNode = nullptr;
    RUNode = nullptr;
    LBNode = nullptr;
    RBNode = nullptr;
    parent = nullptr;
}


// inert point obj 
void pointNode::insertObj(double obj_x, double obj_y, short switch_level, pointAttri attri){
    for (auto field_item: attri.field_string){
        attri_map[field_item.first][field_item.second] += 1;
    }
    node_fcount++;
    
    if (node_level == switch_level + 1){
        Rtree.insert(make_pair(POINT(obj_x, obj_y), attri));
    }

    if (node_level == switch_level + 8){
        return;
    }
    

    // LU sub-region
    if (IfContainPoint(node_x, node_y + node_w/2, node_w/2, obj_x, obj_y)){
        if (!LUNode){
            LUNode = new pointNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        LUNode->insertObj(obj_x, obj_y, switch_level, attri);
        return;
    }
    // RU sub-region
    else if (IfContainPoint(node_x + node_w/2, node_y + node_w/2, node_w/2, obj_x, obj_y)){
        if (!RUNode){
            RUNode = new pointNode(RU, node_x + node_w/2, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        RUNode->insertObj(obj_x, obj_y, switch_level, attri);
        return;
    }
    // LB sub-region
    else if (IfContainPoint(node_x, node_y, node_w/2, obj_x, obj_y)){
        if (!LBNode){
            LBNode = new pointNode(LB, node_x, node_y, node_w/2, node_level + 1, this);
        }
        LBNode->insertObj(obj_x, obj_y, switch_level, attri);
        return;
    }
    // RB sub-reion
    else if (IfContainPoint(node_x + node_w/2, node_y, node_w/2, obj_x, obj_y)){
        if (!RBNode){
            RBNode = new pointNode(RB, node_x + node_w/2, node_y, node_w/2, node_level + 1, this);
        }
        RBNode->insertObj(obj_x, obj_y, switch_level, attri);
        return;
    }
}



pointNode *pointNode::getMinWNode(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level){
    if (node_level == switch_level + 1){
        return this;
    }
    
    // LU sub-region
    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!LUNode){
            return nullptr;
        }
        pointNode *seek = LUNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }

    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RUNode){
            return nullptr;
        }
        pointNode *seek = RUNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!LBNode){
            return nullptr;
        }
        pointNode *seek = LBNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RBNode){
            return nullptr;
        }
        pointNode *seek = RBNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }
    else {
        return this;
    }
}


void pointNode::countField_single(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, map<string, int> &count_map){

    if (node_level == switch_level + 1){
        vector<POINT_ITEM> objs;
        Box box(POINT(box_xMin, box_yMin), POINT(box_xMax, box_yMax));
        Rtree.query(bgi::intersects(box), back_inserter(objs));
        count_map["count"] += objs.size();
        return;
    }

    if (LUNode){
        if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                count_map["count"] += LUNode->node_fcount;
            }
            else{
                LUNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
            }
        }
    }

    if (RUNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                count_map["count"] += RUNode->node_fcount;
            }
            else{
                RUNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
            }
        }
    }

    if (LBNode){
        if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                count_map["count"] += LBNode->node_fcount;
            }
            else{
                LBNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
            }
        }
    }

    if (RBNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                count_map["count"] += RBNode->node_fcount;
            }
            else{
                RBNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
            }
        }
    }
    return;
}


void pointNode::countField_classified(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, map<string, int> &count_map){

    if (node_level == switch_level + 1){
        vector<POINT_ITEM> objs;
        Box box(POINT(box_xMin, box_yMin), POINT(box_xMax, box_yMax));
        Rtree.query(bgi::intersects(box), back_inserter(objs));
        for (auto obj: objs){
            pointAttri field_attri = obj.second;
            string field_item_name = field_attri.field_string[field_name];
            count_map[field_item_name] += 1;
        }
        return;
    }

    if (LUNode){
        if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: LUNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                LUNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
            }
        }
    }

    if (RUNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: RUNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                RUNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
            }
        }
    }

    if (LBNode){
        if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: LBNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                LBNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
            }
        }
    }

    if (RBNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: RBNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                RBNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
            }
        }
    }
    return;
}


void pointNode::countField_ruled(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, vector<string> rules, map<string, int> &count_map){

    if (node_level == switch_level + 1){
        vector<POINT_ITEM> objs;
        Box box(POINT(box_xMin, box_yMin), POINT(box_xMax, box_yMax));
        Rtree.query(bgi::intersects(box), back_inserter(objs));
        for (auto obj: objs){
            pointAttri field_attri = obj.second;
            string field_item_name = field_attri.field_string[field_name];
            map<string, int> field_item = {{field_item_name, 1}};
            updateCount_ruled(field_item, rules, count_map);
        }
        return;
    }

    if (LUNode){
        if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(LUNode->attri_map[field_name], rules, count_map);
            }
            else{
                LUNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
            }
        }
    }

    if (RUNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(RUNode->attri_map[field_name], rules, count_map);
            }
            else{
                RUNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
            }
        }
    }

    if (LBNode){
        if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(LBNode->attri_map[field_name], rules, count_map);
            }
            else{
                LBNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
            }
        }
    }

    if (RBNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(RBNode->attri_map[field_name], rules, count_map);
            }
            else{
                RBNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
            }
        }
    }
    return;
}


void pointNode::getMaxFcount(int z, int &max_fcount){
    if (node_level == (z+8)){
        if (node_fcount > max_fcount){
            max_fcount = node_fcount;
        }
        return;
    }

    if (LUNode){
        LUNode->getMaxFcount(z, max_fcount);
    }
    if (RUNode){
        RUNode->getMaxFcount(z, max_fcount);
    }
    if (LBNode){
        LBNode->getMaxFcount(z, max_fcount);
    }
    if (RBNode){
        RBNode->getMaxFcount(z, max_fcount);
    }
}

void pointNode::getMaxFcount_range(int z, int &max_fcount, double box[]){
    if (node_level == (z+8) && IfOverlapMBB(node_x, node_y, node_w, box[0], box[1], box[2], box[3])){
        if (node_fcount > max_fcount){
            max_fcount = node_fcount;
        }
        return;
    }

    if (LUNode){
        LUNode->getMaxFcount_range(z, max_fcount, box);
    }
    if (RUNode){
        RUNode->getMaxFcount_range(z, max_fcount, box);
    }
    if (LBNode){
        LBNode->getMaxFcount_range(z, max_fcount, box);
    }
    if (RBNode){
        RBNode->getMaxFcount_range(z, max_fcount, box);
    }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------- //
// **********************************************    linestring node functions    *************************************************************** //
// ---------------------------------------------------------------------------------------------------------------------------------------------- //
linestringNode::linestringNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, linestringNode *_parent){
    this->node_type = _node_type;
    this->node_x = _node_x;
    this->node_y = _node_y;
    this->node_w = _node_w;
    this->node_level = _node_level;
    this->node_flen = 0;
    this->parent = _parent;
    this->LUNode = nullptr;
    this->RUNode = nullptr;
    this->LBNode = nullptr;
    this->RBNode = nullptr;
}


linestringNode::~linestringNode(){
    delete LUNode;
    delete RUNode;
    delete LBNode;
    delete RBNode;
    delete parent;
    LUNode = nullptr;
    RUNode = nullptr;
    LBNode = nullptr;
    RBNode = nullptr;
    parent = nullptr;
}


// insert TMBB of a linestring feature
linestringNode *linestringNode::insertTMBB(double obj_xMin, double obj_yMin, double obj_xMax, double obj_yMax, int seg_count, double len, short switch_level, linestringAttri attri){

    for (auto field_item: attri.field_string){
        attri_map[field_item.first][field_item.second] += len;
    }
    node_flen += len;
    
    // LU sub-region
    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= switch_level)){
        if (!LUNode){
            LUNode = new linestringNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        linestringNode *seek = LUNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, seg_count, len, switch_level, attri);
        return seek;
    }

    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= switch_level)){
        if (!RUNode){
            RUNode = new linestringNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level+1, this);
        }
        linestringNode *seek = RUNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, seg_count, len, switch_level, attri);
        return seek;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= switch_level)){
        if (!LBNode){
            LBNode = new linestringNode(LB, node_x, node_y, node_w/2, node_level+1, this);
        }
        linestringNode *seek = LBNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, seg_count, len, switch_level, attri);
        return seek;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= switch_level)){
        if (!RBNode){
            RBNode = new linestringNode(RB, node_x + node_w/2, node_y, node_w/2, node_level+1, this);
        }
        linestringNode *seek = RBNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, seg_count, len, switch_level, attri);
        return seek;
    }
    else {
        return this;
    }
}


// insert the segement of a linestring feature
void linestringNode::insertSeg(SEG seg, short switch_level, linestringAttri attri){
    if (node_level == switch_level + 9){
        return;
    }

    // when the node contain the segement
    // LU sub-region
    Box LU_box = Box(POINT(node_x, node_y+node_w/2), POINT(node_x+node_w/2, node_y+node_w));
    if (bg::intersects(seg, LU_box)){
        if (!LUNode){
            LUNode = new linestringNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        
        bgm::multi_linestring<LINESTRING> intesect_seg;
        bg::intersection(seg, LU_box, intesect_seg);
        double len = bg::length(intesect_seg)/1000;
        for (auto field_item: attri.field_string){
            LUNode->attri_map[field_item.first][field_item.second] += len;
        }
        LUNode->node_flen += len;

        LUNode->insertSeg(seg, switch_level, attri);
    }
    // RU sub-region
    Box RU_box = Box(POINT(node_x+node_w/2, node_y+node_w/2), POINT(node_x+node_w, node_y+node_w));
    if (bg::intersects(seg, RU_box)){
        if (!RUNode){
            RUNode = new linestringNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level + 1, this);
        }

        bgm::multi_linestring<LINESTRING> intesect_seg;
        bg::intersection(seg, RU_box, intesect_seg);
        double len = bg::length(intesect_seg)/1000;
        for (auto field_item: attri.field_string){
            RUNode->attri_map[field_item.first][field_item.second] += len;
        }
        RUNode->node_flen += len;

        RUNode->insertSeg(seg, switch_level, attri);
    }
    // LB sub-region
    Box LB_box = Box(POINT(node_x, node_y), POINT(node_x+node_w/2, node_y+node_w/2));
    if (bg::intersects(seg, LB_box)){
        if (!LBNode){
            LBNode = new linestringNode(LB, node_x, node_y, node_w/2, node_level + 1, this);
        }
        
        bgm::multi_linestring<LINESTRING> intesect_seg;
        bg::intersection(seg, LB_box, intesect_seg);
        double len = bg::length(intesect_seg)/1000;
        for (auto field_item: attri.field_string){
            LBNode->attri_map[field_item.first][field_item.second] += len;
        }
        LBNode->node_flen += len;

        LBNode->insertSeg(seg, switch_level, attri);
    }
    // RB sub-region
    Box RB_box = Box(POINT(node_x+node_w/2, node_y), POINT(node_x+node_w, node_y+node_w/2));
    if (bg::intersects(seg, RB_box)){
        if (!RBNode){
            RBNode = new linestringNode(RB, node_x + node_w/2, node_y, node_w/2, node_level + 1, this);
        }

        bgm::multi_linestring<LINESTRING> intesect_seg;
        bg::intersection(seg, RB_box, intesect_seg);
        double len = bg::length(intesect_seg)/1000;
        for (auto field_item: attri.field_string){
            RBNode->attri_map[field_item.first][field_item.second] += len;
        }
        RBNode->node_flen += len;

        RBNode->insertSeg(seg, switch_level, attri);
    }
    return;
}

void linestringNode::insertFeature(LINESTRING line_g, short switch_level, LINESTRING_ITEM line){
    if (node_level == switch_level + 1){
        Rtree.insert(line);
        return;
    }

    // LU sub-region
    if (LUNode){
        if (bg::intersects(line_g, Box(POINT(node_x, node_y+node_w/2), POINT(node_x+node_w/2, node_y+node_w)))){
            LUNode->insertFeature(line_g, switch_level, line);
        }
    }
    // RU sub-region
    if (RUNode){
        if (bg::intersects(line_g, Box(POINT(node_x+node_w/2, node_y+node_w/2), POINT(node_x+node_w, node_y+node_w)))){
            RUNode->insertFeature(line_g, switch_level, line);
        }
    }
    // LB sub-region
    if (LBNode){
        if (bg::intersects(line_g, Box(POINT(node_x, node_y), POINT(node_x+node_w/2, node_y+node_w/2)))){
            LBNode->insertFeature(line_g, switch_level, line);
        }
    }
    // RB sub-region
    if (RBNode){
        if (bg::intersects(line_g, Box(POINT(node_x+node_w/2, node_y), POINT(node_x+node_w, node_y+node_w/2)))){
            RBNode->insertFeature(line_g, switch_level, line);
        }
    }
    return;
}

linestringNode *linestringNode::getMinWNode(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level){
    if (node_level == switch_level + 1){
        return this;
    }

    // LU sub-region
    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!LUNode){
            return nullptr;
        }
        linestringNode *seek = LUNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }

    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RUNode){
            return nullptr;
        }
        linestringNode *seek = RUNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!LBNode){
            return nullptr;
        }
        linestringNode *seek = LBNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RBNode){
            return nullptr;
        }
        linestringNode *seek = RBNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }
    else {
        return this;
    }
}


void linestringNode::countField_single(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, map<string, double> &count_map){

    if (node_level == switch_level + 1){
        double ibox_xMin = node_x > box_xMin ? node_x : box_xMin;
        double ibox_xMax = (node_x + node_w) > box_xMax ? box_xMax : (node_x + node_w);
        double ibox_yMin = node_y > box_yMin ? node_y : box_yMin;
        double ibox_yMax = (node_y + node_w) > box_yMax ? box_yMax : (node_y + node_w);
        vector<LINESTRING_ITEM> objs;
        Box box(POINT(ibox_xMin, ibox_yMin), POINT(ibox_xMax, ibox_yMax));
        Rtree.query(bgi::intersects(box), back_inserter(objs));
        for (auto obj: objs){
            LINE line_coord = get<1>(obj);
            LINESTRING line;
            for (int i = 0; i < line_coord.size(); i+=2){
                line.push_back(POINT(line_coord[i], line_coord[i+1]));
            }

            bgm::multi_linestring<LINESTRING> intesect_seg;
            bg::intersection(line, box, intesect_seg);
            count_map["length"] += bg::length(intesect_seg)/1000;
        }
        return;
    }
    

    if (LUNode){
        if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                count_map["length"] += LUNode->node_flen;
            }
            else{
                LUNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
            }
        }
    }

    if (RUNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                count_map["length"] += RUNode->node_flen;
            }
            else{
                RUNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
            }
        }
    }

    if (LBNode){
        if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                count_map["length"] += LBNode->node_flen;
            }
            else{
                LBNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
            }
        }
    }

    if (RBNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                count_map["length"] += RBNode->node_flen;
            }
            else{
                RBNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
            }
        }
    }
    return;
}

void linestringNode::countField_classified(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, map<string, double> &count_map){

    if (node_level == switch_level + 1){
        double ibox_xMin = node_x > box_xMin ? node_x : box_xMin;
        double ibox_xMax = (node_x + node_w) > box_xMax ? box_xMax : (node_x + node_w);
        double ibox_yMin = node_y > box_yMin ? node_y : box_yMin;
        double ibox_yMax = (node_y + node_w) > box_yMax ? box_yMax : (node_y + node_w);
        vector<LINESTRING_ITEM> objs;
        Box box(POINT(ibox_xMin, ibox_yMin), POINT(ibox_xMax, ibox_yMax));
        Rtree.query(bgi::intersects(box), back_inserter(objs));
        for (auto obj: objs){
            LINE line_coord = get<1>(obj);
            linestringAttri field_attri = get<2>(obj);
            string field_item_name = field_attri.field_string[field_name];

            LINESTRING line;
            for (int i = 0; i < line_coord.size(); i+=2){
                line.push_back(POINT(line_coord[i], line_coord[i+1]));
            }

            bgm::multi_linestring<LINESTRING> intesect_seg;
            bg::intersection(line, box, intesect_seg);
            count_map[field_item_name] += bg::length(intesect_seg)/1000;
        }
        return;
    }

    if (LUNode){
        if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: LUNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                LUNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
            }
        }
    }

    if (RUNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: RUNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                RUNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
            }
        }
    }

    if (LBNode){
        if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: LBNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                LBNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
            }
        }
    }

    if (RBNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: RBNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                RBNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
            }
        }
    }
    return;
}

void linestringNode::countField_ruled(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, vector<string> rules, map<string, double> &count_map){

    if (node_level == switch_level + 1){
        double ibox_xMin = node_x > box_xMin ? node_x : box_xMin;
        double ibox_xMax = (node_x + node_w) > box_xMax ? box_xMax : (node_x + node_w);
        double ibox_yMin = node_y > box_yMin ? node_y : box_yMin;
        double ibox_yMax = (node_y + node_w) > box_yMax ? box_yMax : (node_y + node_w);
        vector<LINESTRING_ITEM> objs;
        Box box(POINT(ibox_xMin, ibox_yMin), POINT(ibox_xMax, ibox_yMax));
        Rtree.query(bgi::intersects(box), back_inserter(objs));
        for (auto obj: objs){
            LINE line_coord = get<1>(obj);
            linestringAttri field_attri = get<2>(obj);
            string field_item_name = field_attri.field_string[field_name];

            LINESTRING line;
            for (int i = 0; i < line_coord.size(); i+=2){
                line.push_back(POINT(line_coord[i], line_coord[i+1]));
            }

            bgm::multi_linestring<LINESTRING> intesect_seg;
            bg::intersection(line, box, intesect_seg);
            map<string, double> field_item = {{field_item_name, bg::length(intesect_seg)/1000}};
            updateCount_ruled(field_item, rules, count_map);
        }
        return;
    }


    if (LUNode){
        if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(LUNode->attri_map[field_name], rules, count_map);
            }
            else{
                LUNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
            }
        }
    }

    if (RUNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(RUNode->attri_map[field_name], rules, count_map);
            }
            else{
                RUNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
            }
        }
    }

    if (LBNode){
        if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(LBNode->attri_map[field_name], rules, count_map);
            }
            else{
                LBNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
            }
        }
    }

    if (RBNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(RBNode->attri_map[field_name], rules, count_map);
            }
            else{
                RBNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
            }
        }
    }
    return;
}

void linestringNode::getMaxFcount(int z, int &max_fcount){
    if (node_level == z+8){
        if (node_flen > max_fcount){
            max_fcount = node_flen;
        }
        return;
    }

    if (LUNode){
        LUNode->getMaxFcount(z, max_fcount);
    }
    if (RUNode){
        RUNode->getMaxFcount(z, max_fcount);
    }
    if (LBNode){
        LBNode->getMaxFcount(z, max_fcount);
    }
    if (RBNode){
        RBNode->getMaxFcount(z, max_fcount);
    }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------- //
// ***************************************************    polygon node functions    *************************************************************** //
// ---------------------------------------------------------------------------------------------------------------------------------------------- //
polygonNode::polygonNode(nodeType _node_type, double _node_x, double _node_y, double _node_w, unsigned short _node_level, polygonNode *_parent){
    this->node_type = _node_type;
    this->topo_type = 'i';
    this->node_x = _node_x;
    this->node_y = _node_y;
    this->node_w = _node_w;
    this->node_level = _node_level;
    this->node_farea = 0;
    this->parent = _parent;
    this->LUNode = nullptr;
    this->RUNode = nullptr;
    this->LBNode = nullptr;
    this->RBNode = nullptr;
}


polygonNode::~polygonNode(){
    delete LUNode;
    delete RUNode;
    delete LBNode;
    delete RBNode;
    delete parent;
    LUNode = nullptr;
    RUNode = nullptr;
    LBNode = nullptr;
    RBNode = nullptr;
    parent = nullptr;
}


void withinNodeUpdate(polygonNode *node, polygonAttri attri, short switch_level){
    node->topo_type = 'w';
    double area = node->node_w * node->node_w / 1000000;
    for (auto field_item: attri.field_string){
        node->attri_map[field_item.first][field_item.second] += area;
    }
}

void overlapNodeUpdate(polygonNode *node, POLYGON poly, double farea, polygonAttri attri, short switch_level){

    bgm::multi_polygon<POLYGON> intesect_poly;
    Box box(POINT(node->node_x, node->node_y), POINT(node->node_x+node->node_w, node->node_y+node->node_w));
    bg::intersection(poly, box, intesect_poly);
    double i_area = abs(bg::area(intesect_poly) / 1000000);

    if (i_area <= farea){
        for (auto field_item: attri.field_string){
            node->attri_map[field_item.first][field_item.second] += i_area;
        }
        node->node_farea += i_area;
    }
    else{
        double m_area = node->node_w * node->node_w / 1000000 - i_area;
        for (auto field_item: attri.field_string){
            node->attri_map[field_item.first][field_item.second] += m_area;
        }
        node->node_farea += m_area;
    }
}

double calWithinArea(double node_x, double node_y, double node_w, double box_xMin, double box_yMin, double box_xMax, double box_yMax){
    double ibox_xMin = node_x > box_xMin ? node_x : box_xMin;
    double ibox_xMax = (node_x + node_w) > box_xMax ? box_xMax : (node_x + node_w);
    double ibox_yMin = node_y > box_yMin ? node_y : box_yMin;
    double ibox_yMax = (node_y + node_w) > box_yMax ? box_yMax : (node_y + node_w);

    return (ibox_xMax - ibox_xMin) * (ibox_yMax - ibox_yMin) / 1000000;
}

// insert TMBB of a polygon feature
polygonNode * polygonNode::insertTMBB(float obj_xMin, float obj_yMin, float obj_xMax, float obj_yMax, double farea, short switch_level, polygonAttri attri){
    for (auto field_item: attri.field_string){
        attri_map[field_item.first][field_item.second] += farea;
    }
    node_farea += farea;
    
    // LU sub-region
    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= switch_level)){
        if (!LUNode){
            LUNode = new polygonNode(LU, node_x, node_y + node_w/2, node_w/2, node_level+1, this);
        }
        
        polygonNode *seek = LUNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, farea, switch_level, attri);
        return seek;
    }
    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= switch_level)){
        if (!RUNode){
            RUNode = new polygonNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level+1, this);
        }
        
        polygonNode *seek = RUNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, farea, switch_level, attri);
        return seek;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= switch_level)){
        if (!LBNode){
            LBNode = new polygonNode(LB, node_x, node_y, node_w/2, node_level+1, this);
        }
        
        polygonNode *seek = LBNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, farea, switch_level, attri);
        return seek;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, obj_xMin, obj_yMin, obj_xMax, obj_yMax) && (node_level+1 <= switch_level)){
        if (!RBNode){
            RBNode = new polygonNode(RB, node_x+node_w/2, node_y, node_w/2, node_level+1, this);
        }
        
        polygonNode *seek = RBNode->insertTMBB(obj_xMin, obj_yMin, obj_xMax, obj_yMax, farea, switch_level, attri);
        return seek;
    }
    else {
        return this;
    }
}


// insert the edge of a polygon feature
void polygonNode::insertEdge(float edge_xMin, float edge_yMin, float edge_xMax, float edge_yMax, string edge_type, short switch_level, polygonAttri attri){ 
    
    if (node_level == switch_level + 9){
        return;
    }

    // when the node contain the edge
    // LU sub-region

    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax)){
        if (!LUNode){
            LUNode = new polygonNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
        }
        // LUNode->node_fcount++;
        LUNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, switch_level, attri);
        return;
    }
    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax)){
        if (!RUNode){
            RUNode = new polygonNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level+1, this);
        }
        // RUNode->node_fcount++;
        RUNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, switch_level, attri);
        return;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax)){
        if (!LBNode){
            LBNode = new polygonNode(LB, node_x, node_y, node_w/2, node_level+1, this);
        }
        // LBNode->node_fcount++;
        LBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, switch_level, attri);
        return;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax)){
        if (!RBNode){
            RBNode = new polygonNode(RB, node_x+node_w/2, node_y, node_w/2, node_level+1, this);
        }
        // RBNode->node_fcount++;
        RBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, switch_level, attri);
        return;
    }
    // when the node overlap with the edge
    else{
        if (IfInterctSeg(node_x, node_y + node_w/2, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type)){
            if (!LUNode){
                LUNode = new polygonNode(LU, node_x, node_y + node_w/2, node_w/2, node_level + 1, this);
            }
            // LUNode->node_fcount++;
            LUNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, switch_level, attri);
        }
        if (IfInterctSeg(node_x + node_w/2, node_y+node_w/2, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type)){
            if (!RUNode){
                RUNode = new polygonNode(RU, node_x + node_w/2, node_y+node_w/2, node_w/2, node_level + 1, this);
            }
            // RUNode->node_fcount++;
            RUNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, switch_level, attri);
        }
        if (IfInterctSeg(node_x, node_y, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type)){
            if (!LBNode){
                LBNode = new polygonNode(LB, node_x, node_y, node_w/2, node_level + 1, this);
            }
            // LBNode->node_fcount++;
            LBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, switch_level, attri);
        }
        if (IfInterctSeg(node_x + node_w/2, node_y, node_w/2, edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type)){
            if (!RBNode){
                RBNode = new polygonNode(RB, node_x + node_w/2, node_y, node_w/2, node_level + 1, this);
            }
            // RBNode->node_fcount++;
            RBNode->insertEdge(edge_xMin, edge_yMin, edge_xMax, edge_yMax, edge_type, switch_level, attri);
        }
    }
    return;
}

// insert the polygon feature
void polygonNode::insertFeature(float box_xMin, float box_yMin, float box_xMax, float box_yMax, short switch_level, float ring_x[], float ring_y[], int len, polygonAttri attri, POLYGON poly, double farea, POLYGON_ITEM ring){
    if (topo_type == 'i'){
        if (node_level == switch_level + 1){
            Rtree.insert(ring);
        }
    }
    else if (topo_type == 'w'){
        return;
    }


    if (node_level == switch_level + 9){
        return;
    }



    // LU sub-region
    if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        // if node is not exist, the node must within the polygon feature
        if (!LUNode){
            if (IfwithinBox(node_x+node_w/4, node_y+3*node_w/4, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (IfwithinPolygon(node_x+node_w/4, node_y+3*node_w/4, ring_x, ring_y, len)){
                    LUNode = new polygonNode(LU, node_x, node_y+node_w/2, node_w/2, node_level+1, this);
                    withinNodeUpdate(LUNode, attri, switch_level);
                }
            }
        }
        // if node exists
        else{
            if (LUNode->topo_type == 'i'){
                overlapNodeUpdate(LUNode, poly, farea, attri, switch_level);
                LUNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, ring_x, ring_y, len, attri, poly, farea, ring);
            }
        }
    }
    // RU sub-region
    if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RUNode){
            if (IfwithinBox(node_x+3*node_w/4, node_y+3*node_w/4, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (IfwithinPolygon(node_x+3*node_w/4, node_y+3*node_w/4, ring_x, ring_y, len)){
                    RUNode = new polygonNode(RU, node_x+node_w/2, node_y+node_w/2, node_w/2, node_level+1, this);
                    withinNodeUpdate(RUNode, attri, switch_level);
                }
            }
        }
        else{
            if (RUNode->topo_type == 'i'){
                overlapNodeUpdate(RUNode, poly, farea, attri, switch_level);
                RUNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, ring_x, ring_y, len, attri, poly, farea, ring);
            }
        } 
    }
    // LB sub-region
    if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!LBNode){
            if (IfwithinBox(node_x+node_w/4, node_y+node_w/4, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (IfwithinPolygon(node_x+node_w/4, node_y+node_w/4, ring_x, ring_y, len)){
                    LBNode = new polygonNode(LB, node_x, node_y, node_w/2, node_level+1, this);
                    withinNodeUpdate(LBNode, attri, switch_level);
                }  
            }
        }
        else{
            if (LBNode->topo_type == 'i'){
                overlapNodeUpdate(LBNode, poly, farea, attri, switch_level);
                LBNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, ring_x, ring_y, len, attri, poly, farea, ring);
            }
        }
    }
    // RB sub-region
    if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RBNode){
            if (IfwithinBox(node_x+3*node_w/4, node_y+node_w/4, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (IfwithinPolygon(node_x+3*node_w/4, node_y+node_w/4, ring_x, ring_y, len)){
                    RBNode = new polygonNode(RB, node_x+node_w/2, node_y, node_w/2, node_level+1, this);
                    withinNodeUpdate(RBNode, attri, switch_level);
                } 
            }
        }
        else{
            if (RBNode->topo_type == 'i'){
                overlapNodeUpdate(RBNode, poly, farea, attri, switch_level);
                RBNode->insertFeature(box_xMin, box_yMin, box_xMax, box_yMax,switch_level, ring_x, ring_y, len, attri, poly, farea, ring);
            }
        }
    }
    return; 
}


polygonNode *polygonNode::getMinWNode(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level){
    if (node_level == switch_level + 1){
        return this;
    }

    // LU sub-region
    if (IfContainMBB(node_x, node_y + node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!LUNode){
            return nullptr;
        }
        polygonNode *seek = LUNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }

    // RU sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RUNode){
            return nullptr;
        }
        polygonNode *seek = RUNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }
    // LB sub-region
    else if (IfContainMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!LBNode){
            return nullptr;
        }
        polygonNode *seek = LBNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }
    // RB sub-region
    else if (IfContainMBB(node_x + node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
        if (!RBNode){
            return nullptr;
        }
        polygonNode *seek = RBNode->getMinWNode(box_xMin, box_yMin, box_xMax, box_yMax, switch_level);
        return seek;
    }
    else {
        return this;
    }
}



void polygonNode::countField_single(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, map<string, double> &count_map){

    if (node_level == switch_level + 1){
        double ibox_xMin = node_x > box_xMin ? node_x : box_xMin;
        double ibox_xMax = (node_x + node_w) > box_xMax ? box_xMax : (node_x + node_w);
        double ibox_yMin = node_y > box_yMin ? node_y : box_yMin;
        double ibox_yMax = (node_y + node_w) > box_yMax ? box_yMax : (node_y + node_w);

        vector<POLYGON_ITEM> objs;
        Box box(POINT(ibox_xMin, ibox_yMin), POINT(ibox_xMax, ibox_yMax));
        Rtree.query(bgi::intersects(box), back_inserter(objs));
        for (auto obj: objs){
            RING ring_coord = get<1>(obj);
            POLYGON poly;
            for (int i = 0; i < ring_coord.size(); i+=2){
                poly.outer().push_back(POINT(ring_coord[i], ring_coord[i+1]));
            }

            bgm::multi_polygon<POLYGON> intesect_poly;
            bg::intersection(poly, box, intesect_poly);
            count_map["area"] += abs(bg::area(intesect_poly) / 1000000);
        }
        return;
    }
    

    if (LUNode){
        if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (LUNode->topo_type == 'i'){
                    count_map["area"] += LUNode->node_farea;
                }
                else{
                    count_map["area"] += LUNode->node_w * LUNode->node_w / 1000000;
                }
            }
            else{
                if (LUNode->topo_type == 'i'){
                    LUNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
                }
                else{
                    count_map["area"] += calWithinArea(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                }
            }
        }
    }

    if (RUNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (RUNode->topo_type == 'i'){
                    count_map["area"] += RUNode->node_farea;
                }
                else{
                    count_map["area"] += RUNode->node_w * RUNode->node_w / 1000000;
                }
            }
            else{
                if (RUNode->topo_type == 'i'){
                    RUNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
                }
                else{
                    count_map["area"] += calWithinArea(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                }
            }
        }
    }

    if (LBNode){
        if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (LBNode->topo_type == 'i'){
                    count_map["area"] += LBNode->node_farea;
                }
                else{
                    count_map["area"] += LBNode->node_w * LBNode->node_w / 1000000;
                }
            }
            else{
                if (LBNode->topo_type == 'i'){
                    LBNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
                }
                else{
                    count_map["area"] += calWithinArea(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                }
            }
        }
    }

    if (RBNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                if (RBNode->topo_type == 'i'){
                    count_map["area"] += RBNode->node_farea;
                }
                else{
                    count_map["area"] += RBNode->node_w * RBNode->node_w / 1000000;
                }
            }
            else{
                if (RBNode->topo_type == 'i'){
                    RBNode->countField_single(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, count_map);
                }
                else{
                    count_map["area"] += calWithinArea(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                }
            }
        }
    }
    return;
}


void polygonNode::countField_classified(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, map<string, double> &count_map){
    if (node_level == switch_level + 1){
        double ibox_xMin = node_x > box_xMin ? node_x : box_xMin;
        double ibox_xMax = (node_x + node_w) > box_xMax ? box_xMax : (node_x + node_w);
        double ibox_yMin = node_y > box_yMin ? node_y : box_yMin;
        double ibox_yMax = (node_y + node_w) > box_yMax ? box_yMax : (node_y + node_w);

        vector<POLYGON_ITEM> objs;
        Box box(POINT(ibox_xMin, ibox_yMin), POINT(ibox_xMax, ibox_yMax));
        Rtree.query(bgi::intersects(box), back_inserter(objs));
        for (auto obj: objs){
            RING ring_coord = get<1>(obj);
            polygonAttri field_attri = get<2>(obj);
            string field_item_name = field_attri.field_string[field_name];
            
            POLYGON poly;
            for (int i = 0; i < ring_coord.size(); i+=2){
                poly.outer().push_back(POINT(ring_coord[i], ring_coord[i+1]));
            }

            bgm::multi_polygon<POLYGON> intesect_poly;
            bg::intersection(poly, box, intesect_poly);
            count_map[field_item_name] += abs(bg::area(intesect_poly) / 1000000);
        }
        return;
    }
    

    if (LUNode){
        if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: LUNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                if (LUNode->topo_type == 'i'){
                    LUNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
                }
                else{
                    for (auto field_item: LUNode->attri_map[field_name]){
                        count_map[field_item.first] += calWithinArea(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                    }
                }
            }
        }
    }

    if (RUNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: RUNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                if (RUNode->topo_type == 'i'){
                    RUNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
                }
                else{
                    for (auto field_item: RUNode->attri_map[field_name]){
                        count_map[field_item.first] += calWithinArea(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                    }
                }
            }
        }
    }

    if (LBNode){
        if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: LBNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                if (LBNode->topo_type == 'i'){
                    LBNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
                }
                else{
                    for (auto field_item: LBNode->attri_map[field_name]){
                        count_map[field_item.first] += calWithinArea(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                    }
                }
            }
        }
    }

    if (RBNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                for (auto field_item: RBNode->attri_map[field_name]){
                    count_map[field_item.first] += field_item.second;
                }
            }
            else{
                if (RBNode->topo_type == 'i'){
                    RBNode->countField_classified(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, count_map);
                }
                else{
                    for (auto field_item: RBNode->attri_map[field_name]){
                        count_map[field_item.first] += calWithinArea(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                    }
                }
            }
        }
    }
    return;
}


void polygonNode::countField_ruled(double box_xMin, double box_yMin, double box_xMax, double box_yMax, short switch_level, string field_name, vector<string> rules, map<string, double> &count_map){
    if (node_level == switch_level + 1){
        double ibox_xMin = node_x > box_xMin ? node_x : box_xMin;
        double ibox_xMax = (node_x + node_w) > box_xMax ? box_xMax : (node_x + node_w);
        double ibox_yMin = node_y > box_yMin ? node_y : box_yMin;
        double ibox_yMax = (node_y + node_w) > box_yMax ? box_yMax : (node_y + node_w);

        vector<POLYGON_ITEM> objs;
        Box box(POINT(ibox_xMin, ibox_yMin), POINT(ibox_xMax, ibox_yMax));
        Rtree.query(bgi::intersects(box), back_inserter(objs));
        for (auto obj: objs){
            RING ring_coord = get<1>(obj);
            polygonAttri field_attri = get<2>(obj);
            string field_item_name = field_attri.field_string[field_name];
            
            POLYGON poly;
            for (int i = 0; i < ring_coord.size(); i+=2){
                poly.outer().push_back(POINT(ring_coord[i], ring_coord[i+1]));
            }

            bgm::multi_polygon<POLYGON> intesect_poly;
            bg::intersection(poly, box, intesect_poly);
            
            map<string, double> field_item = {{field_item_name, abs(bg::area(intesect_poly)/1000000)}};
            updateCount_ruled(field_item, rules, count_map);
        }
        return;
    }

    if (LUNode){
        if (IfOverlapMBB(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(LUNode->attri_map[field_name], rules, count_map);
            }
            else{
                if (LUNode->topo_type == 'i'){
                    LUNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
                }
                else{
                    double i_area = calWithinArea(node_x, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                    updateCount_withinNode_ruled(LUNode->attri_map[field_name], i_area, rules, count_map);
                }
            }
        }
    }

    if (RUNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(RUNode->attri_map[field_name], rules, count_map);
            }
            else{
                if (RUNode->topo_type == 'i'){
                    RUNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
                }
                else{
                    double i_area = calWithinArea(node_x+node_w/2, node_y+node_w/2, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                    updateCount_withinNode_ruled(RUNode->attri_map[field_name], i_area, rules, count_map);
                }
            }
        }
    }

    if (LBNode){
        if (IfOverlapMBB(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(LBNode->attri_map[field_name], rules, count_map);
            }
            else{
                if (LBNode->topo_type == 'i'){
                    LBNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
                }
                else{
                    double i_area = calWithinArea(node_x, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                    updateCount_withinNode_ruled(LBNode->attri_map[field_name], i_area, rules, count_map);
                }
            }
        }
    }

    if (RBNode){
        if (IfOverlapMBB(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
            if (IfwithinBox(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax)){
                updateCount_ruled(RBNode->attri_map[field_name], rules, count_map);
            }
            else{
                if (RBNode->topo_type == 'i'){
                    RBNode->countField_ruled(box_xMin, box_yMin, box_xMax, box_yMax, switch_level, field_name, rules, count_map);
                }
                else{
                    double i_area = calWithinArea(node_x+node_w/2, node_y, node_w/2, box_xMin, box_yMin, box_xMax, box_yMax);
                    updateCount_withinNode_ruled(RBNode->attri_map[field_name], i_area, rules, count_map);
                }
            }
        }
    }
    return;
}