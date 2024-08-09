cmake_minimum_required (VERSION 3.9)

project (HiVQ
  VERSION 1.0.0
  LANGUAGES CXX
)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_ROOT}" ${CMAKE_MODULE_PATH})


find_package(PNG REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)
find_library(MAPNIK_LIB mapnik)
find_library(SHP_LIB shp)


if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS}${OpenMP_EXE_LINKER_FLAGS}")
endif()
set(CMAKE_CXX_FLAGS "-lmapnik -lshp -lmpich")
set(INCLUDE_DIR /usr/local/include /usr/include)
set(CORE_LIB   ${PROJECT_SOURCE_DIR}/bin/libhigis.a  ${PROJECT_SOURCE_DIR}/bin/libhigis.so)


add_definitions(-std=c++17 -Wall -O3 -g -DBIGINT)

include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${INCLUDE_DIR}
                    ${MPI_INCLUDE_PATH}
                    )

add_executable(visualServer  
              ${PROJECT_SOURCE_DIR}/src/buildIndex.cpp
              ${PROJECT_SOURCE_DIR}/src/treeNode.cpp 
              ${PROJECT_SOURCE_DIR}/src/intersectJudge.cpp 
              ${PROJECT_SOURCE_DIR}/src/pointPlot.cpp
              ${PROJECT_SOURCE_DIR}/src/linestringPlot.cpp
              ${PROJECT_SOURCE_DIR}/src/polygonPlot.cpp
              ${PROJECT_SOURCE_DIR}/src/visualServer.cpp
              )
target_link_libraries(visualServer 
${MAPNIK_LIB} ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${CORE_LIB} ${SHP_LIB} ${Python3_LIBRARIES} -lpthread)








