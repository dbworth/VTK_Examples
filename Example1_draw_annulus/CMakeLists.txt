cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(myproject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE RelWithDebInfo) # Debug, Release, RelWithDebInfo

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(draw_annulus draw_annulus.cpp)
target_link_libraries(draw_annulus ${VTK_LIBRARIES})
