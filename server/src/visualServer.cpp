#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "crow.h"
#include "png.h"
#include <regex>
#include "treeNode.hpp"
#include "buildIndex.hpp"
#include "pointPlot.hpp"
#include "linestringPlot.hpp"
#include "polygonPlot.hpp"
#include <shapefil.h>
using namespace std;
#define TILE_SIZE 256
#define L	20037508.34


struct ServerMiddleware
{
	std::string message;

	ServerMiddleware()
	{
		message = "foo";
	}


	void setMessage( std::string newMsg )
	{
		message = newMsg;
	}


	struct context
	{
	};

	void before_handle( crow::request & /*req*/, crow::response & /*res*/, context & /*ctx*/ )
	{
		CROW_LOG_DEBUG << " - MESSAGE: " << message;
	}


	void after_handle( crow::request & /*req*/, crow::response & /*res*/, context & /*ctx*/ )
	{
		/* no-op */
	}
};



bool checkInput(string data_type, string data_path, string data_srs){
    bool dataPath_statue = false;
    bool dataSrs_statue = false;
    if (data_type == "shape"){
        if(data_path.find(".shp") == string::npos){
            cout << "【ERROR】 please input correct file format!" << endl;
            cout << "【ERROR】 " << data_path << endl;
            dataPath_statue = true;
        }
    }
    else if (data_type == "csv"){
        if(data_path.find(".csv") == string::npos){
            cout << "【ERROR】 please input correct file format!" << endl;
            cout << "【ERROR】 " << data_path << endl;
            dataPath_statue = true;
        }
    }
    
    
    if (data_srs != "4326" and data_srs != "3857"){
        cout << "【ERROR】 please input correct srs!" << endl;
        cout << "【ERROR】 " << data_path << " " << data_srs << endl;
        dataSrs_statue = true;
    }
    return (dataPath_statue or dataSrs_statue); 
}


void splitField(string str, char split, FieldSet& res){
    if (str == "null"){
        return;
    }


    istringstream iss(str);
    string token;
    while (getline(iss, token, split)){
        istringstream iss2(token);
        string token2;
        vector<string> infos;
        while (getline(iss2, token2, '-')){
            infos.push_back(token2);
        }
        
        if (infos[1] == "s"){
            res.fields_classify.push_back(infos[0]);
        }
        else if (infos[1] == "i" or infos[1] == "d"){
            res.fields_continuous.push_back(infos[0]);
        }
    }
}

vector<string> mergeFields(FieldSet fields){
    vector<string> all_fields;
    vector<string> fields_s = fields.fields_classify;
    vector<string> fields_r = fields.fields_continuous;
    all_fields.insert(all_fields.end(), fields_s.begin(), fields_s.end());
    all_fields.insert(all_fields.end(), fields_r.begin(), fields_r.end());
    return all_fields;
}



int main(int argc, char **argv){    
    // get parameters
    FieldSet point_fields, linestring_fields, polygon_fields;
    string point_fileType = argv[1];
    string point_shpPath = argv[2];
    string point_srs = argv[3];
    string point_fields_list = argv[4];
    
    string linestring_fileType = argv[5];
    string linestring_shpPath = argv[6];
    string linestring_srs = argv[7];
    string linestring_fields_list = argv[8];
    
    string polygon_fileType = argv[9];
    string polygon_shpPath = argv[10];
    string polygon_srs = argv[11];
    string polygon_fields_list = argv[12];
    
    crow::App<ServerMiddleware> app;
    splitField(point_fields_list, ',', point_fields);
    splitField(linestring_fields_list, ',', linestring_fields);
    splitField(polygon_fields_list, ',', polygon_fields);
    

    // ********************************************************************************************************************* //
	// check input information
	// ********************************************************************************************************************* //
    bool point_statue = checkInput(point_fileType, point_shpPath, point_srs);
    bool linestring_statue = checkInput(linestring_fileType, linestring_shpPath, linestring_srs);
    bool polygon_statue = checkInput(polygon_fileType, polygon_shpPath, polygon_srs);
    if (point_statue or linestring_statue or polygon_statue){
        return 0;
    }
    
    
    

    // ********************************************************************************************************************* //
	// build index
	// ********************************************************************************************************************* //
    cout << "--------------------------------------------------------------------------------" << endl;
    pointNode *pointRNode;
    short point_switchLevel;
    pair<pointNode*, short> point_info_shp;

    if (point_fileType == "shape"){
        point_info_shp = pointIndex_shp(point_shpPath, point_srs, point_fields);
        pointRNode = point_info_shp.first;
        point_switchLevel = point_info_shp.second;
    }
    
    cout << "<<<<<<< read shapefile: " << point_shpPath << " successfully >>>>>>>" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;


    cout << "--------------------------------------------------------------------------------" << endl;
    linestringNode *linestringRNode;
    short linestring_switchLevel;
    pair<linestringNode*, short> linestring_info_shp;

    if (linestring_fileType == "shape"){
        linestring_info_shp = linestringIndex_shp(linestring_shpPath, linestring_srs, linestring_fields);
        linestringRNode = linestring_info_shp.first;
        linestring_switchLevel = linestring_info_shp.second;
    }
    
    cout << "<<<<<<< read shapefile: " << linestring_shpPath << " successfully >>>>>>>" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;


    cout << "--------------------------------------------------------------------------------" << endl;
    polygonNode *polygonRNode;
    short polygon_switchLevel;
    pair<polygonNode*, short> polygon_info_shp;

    if (polygon_fileType == "shape"){
        polygon_info_shp = polygonIndex_shp(polygon_shpPath, polygon_srs, polygon_fields);
        polygonRNode = polygon_info_shp.first;
        polygon_switchLevel = polygon_info_shp.second;
    }
    
    cout << "<<<<<<< read shapefile: " << polygon_shpPath << " successfully >>>>>>>" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;

    
    // ********************************************************************************************************************* //
	// http service of showing meta data of the dataset
	// ********************************************************************************************************************* //
    CROW_ROUTE(app, "/HiVQ/<string>" ).name("getMetaData")
		([&] (const crow::request & req, crow::response & res, string shpType) {
            string data_path, file_type;
            string meta_info = "";
            if (shpType == "point"){
                data_path = point_shpPath;
                file_type = point_fileType;
                meta_info += ("Point/" + point_srs + "/");
            }
            else if (shpType == "line"){
                data_path = linestring_shpPath;
                file_type = linestring_fileType;
                meta_info += ("Linestring/" + linestring_srs + "/");
            }
            else if (shpType == "polygon"){
                data_path = polygon_shpPath;
                file_type = polygon_fileType;
                meta_info += ("Polygon/" + polygon_srs + "/");
            }
            
            std::stringstream ss, ss2, ss3, ss4;
            // file type
            if (file_type == "shape"){
                Shapefile shp;
                GeoData data = shp.Read(data_path);
                
                float xmin = data.extent.xmin;
                float xmax = data.extent.xmax;
                float ymin = data.extent.ymin;
                float ymax = data.extent.ymax;
                ss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << xmin;
                string xmin_str = ss.str();
                ss2 << std::setiosflags(std::ios::fixed) << std::setprecision(1) << ymin;
                string ymin_str = ss2.str();
                ss3 << std::setiosflags(std::ios::fixed) << std::setprecision(1) << xmax;
                string xmax_str = ss3.str();
                ss4 << std::setiosflags(std::ios::fixed) << std::setprecision(1) << ymax;
                string ymax_str = ss4.str();

                if (shpType == "point"){
                    meta_info += (to_string(data.feature_count) + "/");
                }
                else if (shpType == "line"){
                    meta_info += (to_string(data.feature_count) + "," + to_string(linestringRNode->node_flen) + " KM/");
                }
                else if (shpType == "polygon"){
                    meta_info += (to_string(data.feature_count) + "," + to_string(polygonRNode->node_farea) + "KM2/");
                }
                
                meta_info += (xmin_str + "," + ymin_str + "-" + xmax_str + "," + ymax_str + "/");

                DBFHandle data_dbf = DBFOpen(shp.AuxPath(data_path, ".dbf").c_str(), "r");
                char *field = new char[12];
                int *pnWidth, *pnDecimals; 
                string field_type[3] = {"string", "int", "double"};
                for (int j = 0; j < DBFGetFieldCount(data_dbf); j++){
                    auto field_info = DBFGetFieldInfo(data_dbf, j, field, pnWidth, pnDecimals);
                    string field_str = field;
                    meta_info += (field_str + "-" + field_type[field_info] + ", ");
                }
            }

            res.write(meta_info);
            res.set_header("Content-Type", "text/html");
			res.set_header("Access-Control-Allow-Origin", "*");
			res.end();
        });


    // ********************************************************************************************************************* //
	// http service of getting the item of the field
	// ********************************************************************************************************************* //
    CROW_ROUTE(app, "/HiVQ/<string>/<string>/<string>").name("getFieldItem")
        ([&] (const crow::request & req, crow::response & res, string shpType, string field_name, string field_type) {
            string field_info = "";
            if (shpType == "point"){
                vector<string> fields_list = mergeFields(point_fields);
                for (int i = 0; i < fields_list.size(); i++){
                    if (field_name == fields_list[i]){
                        map<string, int> fields = (pointRNode->attri_map)[field_name];
                        for (auto field: fields){
                            field_info += (field.first + "," + to_string(field.second) + ",");
                        }
                    }
                }
            }
            else if (shpType == "line"){
                vector<string> fields_list = mergeFields(linestring_fields);
                for (int i = 0; i < fields_list.size(); i++){
                    if (field_name == fields_list[i]){
                        map<string, double> fields = (linestringRNode->attri_map)[field_name];
                        for (auto field: fields){
                            field_info += (field.first + "," + to_string(field.second) + ",");
                        }
                    }
                }
            }
            else if (shpType == "polygon"){
                vector<string> fields_list = mergeFields(polygon_fields);
                for (int i = 0; i < fields_list.size(); i++){
                    if (field_name == fields_list[i]){
                        map<string, double> fields = (polygonRNode->attri_map)[field_name];
                        for (auto field: fields){
                            field_info += (field.first + "," + to_string(field.second) + ",");
                        }
                    }
                }
            }

            res.write(field_info);
            res.set_header("Content-Type", "text/html");
			res.set_header("Access-Control-Allow-Origin", "*");
			res.end();
        });




    // ********************************************************************************************************************* //
	// WTMS of vector data visualization with attribution 
	// ********************************************************************************************************************* //
    CROW_ROUTE(app, "/HiVQ/<string>/<string>/<int>/<int>/<int>.png").name("vision")
		([&] (const crow::request & req, crow::response & res, string shpType, string style_describer, int z, int x, int y) {
			string img_str;
            // style_describer = "single-255,0,0,1,2";
            // style_describer = "ruled-code-int-=1001,255,0,0,1,2-!=1001,0,255,0,1,2";
            // style_describer = "classified-code-int-1001,255,0,0,1,2-1002,255,0,0,1,2-1003,0,255,0,1,2";
            // style_describer = "heat-jet-5-10,15,20,25";
            try{
                if (shpType == "point"){
                    if (z <= point_switchLevel){
					    img_str = disDrivenPlot(z, x, y, style_describer, pointRNode);
                    }
                    else{
                        img_str = dataDrivenPlot(z, x, y, style_describer, pointRNode, point_switchLevel);
                    }
                }
                else if (shpType == "line"){
                    if (z <= linestring_switchLevel){
                        img_str = disDrivenPlot(z, x, y, style_describer, linestringRNode);
                    }
                    else{
                        img_str = dataDrivenPlot(z, x, y, style_describer, linestringRNode, linestring_switchLevel);
                    }
                }
                else if (shpType == "polygon"){
                    if (z <= polygon_switchLevel){
                        img_str = disDrivenPlot(z, x, y, style_describer, polygonRNode);
                    }
                    else{
                        img_str = dataDrivenPlot(z, x, y, style_describer, polygonRNode, polygon_switchLevel);
                    }
                }

				res.write(img_str);
				res.set_header("Content-Type", "image/png" );
				res.set_header("Access-Control-Allow-Origin", "*");
				res.end();
            }
            catch(const char* msg){
				res.code = 500;
                ostringstream os;
				os << msg << "\n";
				res.write(os.str());
				res.set_header("Content-Type", "text/html");
				res.end();
            }
        });




    // ********************************************************************************************************************* //
	// WTMS of vector data visualization by spatial range search
	// ********************************************************************************************************************* //
    CROW_ROUTE(app, "/HiVQ/<string>/<string>/<double>/<double>/<double>/<double>/<int>/<int>/<int>.png").name("range_vision")
		([&] (const crow::request & req, crow::response & res, string shpType, string style_describer, double xMin, double yMin, double xMax, double yMax, int z, int x, int y) {
			string img_str;
            double range_box[] = {xMin, yMin, xMax, yMax};
            // style_describer = "single-255,0,0,1,2";
            // style_describer = "ruled-code-int-=1001,255,0,0,1,2-!=1001,0,255,0,1,2";
            // style_describer = "classified-code-int-1001,255,0,0,1,2-1002,255,0,0,1,2-1003,0,255,0,1,2";
            // style_describer = "heat-jet-5-10,15,20,25";
            try{
                if (shpType == "point"){
                    if (z <= point_switchLevel){
					    img_str = disDrivenPlot_range(z, x, y, range_box, style_describer, pointRNode);
                    }
                    else{
                        img_str = dataDrivenPlot_range(z, x, y, range_box, style_describer, pointRNode, point_switchLevel);
                    }
                }
                else if (shpType == "line"){
                    if (z <= linestring_switchLevel){
                        img_str = disDrivenPlot_range(z, x, y, range_box, style_describer, linestringRNode);
                    }
                    else{
                        img_str = dataDrivenPlot_range(z, x, y, range_box, style_describer, linestringRNode, linestring_switchLevel);
                    }
                }
                else if (shpType == "polygon"){
                    if (z <= polygon_switchLevel){
                        img_str = disDrivenPlot_range(z, x, y, range_box, style_describer, polygonRNode);
                    }
                    else{
                        img_str = dataDrivenPlot_range(z, x, y, range_box, style_describer, polygonRNode, polygon_switchLevel);
                    }
                }

				res.write(img_str);
				res.set_header("Content-Type", "image/png" );
				res.set_header("Access-Control-Allow-Origin", "*");
				res.end();
            }
            catch(const char* msg){
				res.code = 500;
                ostringstream os;
				os << msg << "\n";
				res.write(os.str());
				res.set_header("Content-Type", "text/html");
				res.end();
            }
        });



    
    // ********************************************************************************************************************* //
	// service of statistic results of spatial range search
	// ********************************************************************************************************************* //
    CROW_ROUTE(app, "/HiVQ/rangeStatistic/<string>/<string>/<double>/<double>/<double>/<double>").name("rangeStatistic")
        ([&] (const crow::request & req, crow::response & res, string shpType, string style_describer, double xMin, double yMin, double xMax, double yMax) {
            // style_describer = "single-255,0,0,1,2";
            // style_describer = "ruled-code-int-=1001,255,0,0,1,2-!=1001,0,255,0,1,2";
            // style_describer = "classified-code-int-1001,255,0,0,1,2-1002,255,0,0,1,2-1003,0,255,0,1,2";

            
            string statistic_str = "";
            
            struct timeval	t1, t2;
            gettimeofday(&t1, NULL);

            if (shpType == "point"){
                statistic_str = fieldStatistic(xMin, yMin, xMax, yMax, style_describer, point_switchLevel, pointRNode);
            }
            else if (shpType == "line"){
                statistic_str = fieldStatistic(xMin, yMin, xMax, yMax, style_describer, linestring_switchLevel, linestringRNode);
            }
            else if (shpType == "polygon"){
                statistic_str = fieldStatistic(xMin, yMin, xMax, yMax, style_describer, polygon_switchLevel, polygonRNode);
            }

            gettimeofday(&t2, NULL);
            float time_use = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000.0;
            cout << "result: " << statistic_str << " , time: " << time_use << " ms" << endl;

            res.write(statistic_str);
            res.set_header("Content-Type", "text/html");
			res.set_header("Access-Control-Allow-Origin", "*");
			res.end();
        });

    
    app.port(10085).multithreaded().run();
}