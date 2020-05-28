cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 14)

project(HyperMandelBulb LANGUAGES C CXX)

find_package(ospray 2.1.1 REQUIRED)
find_package(ospcommon 1.3.1 REQUIRED)
find_package(TBB REQUIRED tbb)
find_package(Open3D 0.9.0 REQUIRED)
find_package(nlohmann_json REQUIRED)

link_directories(${Open3D_LIBRARY_DIRS})

file(GLOB SRCS
    "src/*.h"
    "src/*.cpp"
)

add_executable(HyperMandelBulb ${SRCS})

target_link_libraries(HyperMandelBulb ospray::ospray TBB::tbb nlohmann_json::nlohmann_json ospcommon::ospcommon ${Open3D_LIBRARIES})
target_include_directories(HyperMandelBulb PUBLIC ${Open3D_INCLUDE_DIRS})

#nlohmann_json::nlohmann_json 
#target_include_directories(main PUBLIC ${Open3D_INCLUDE_DIRS})
#target_link_libraries (main ${Open3D_LIBRARIES})
#link_directories(${Open3D_LIBRARY_DIRS})




