# !/bin/bash
#point
point_type="shape"
point_shpPath="./data/China_mainland_poi/China_mainland_poi.shp"
point_srs="4326"
point_field="fclass-s"
#linestring
linestring_type="shape"
linestring_shpPath="./data/China_mainland_railway/China_mainland_railway.shp"
linestring_srs="4326"
linestring_field="null"
#polygon
polygon_type="shape"
polygon_shpPath="./data/China_mainland_traffic/China_mainland_traffic.shp"
polygon_srs="4326"
polygon_field="null"


nohup ./build/visualServer $point_type $point_shpPath $point_srs $point_field $linestring_type $linestring_shpPath $linestring_srs $linestring_field $polygon_type $polygon_shpPath $polygon_srs $polygon_field > ./TMS_server.log 2>&1 &
python3 ../browser/release.py