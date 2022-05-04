cmake_minimum_required(VERSION 3.9)
set(CMAKE_C_COMPILER /opt/intel/oneapi/compiler/2021.2.0/mac/bin/intel64/icc)
set(CMAKE_CXX_COMPILER /opt/intel/oneapi/compiler/2021.2.0/mac/bin/intel64/icpc)
set(CMAKE_CXX_STANDARD 17)

project(MixtureCP)

include_directories(
${CMAKE_SOURCE_DIR}
/Users/nguyenly/local/eigen/3.4.1/eigen-3.4-rc1
/usr/local/Cellar/pugixml/1.10/include
)

link_directories(
/usr/local/Cellar/pugixml/1.10/lib
)

add_executable(
MixtureCP 
main.cpp
Common/commonMacros.h
Math/poly34.cpp Math/poly34.h 
Physics/pengrobinson.cpp Physics/pengrobinson.h 
)


target_link_libraries(
MixtureCP
pugixml
)

target_compile_definitions(
MixtureCP PUBLIC 
NASAPOLYPATH="${CMAKE_SOURCE_DIR}/NASAPOLY/nasa_poly.xml"
)


