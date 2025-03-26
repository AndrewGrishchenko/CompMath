cmake_minimum_required(VERSION 3.10.0)
project(cm1 VERSION 0.1.0)

find_package(Eigen3 REQUIRED)
include_directories(${Eigen3_INCLUDE_DIRS})

file(GLOB SOURCES *.cpp)

add_executable(cm1 ${SOURCES})
target_link_libraries(cm1 Eigen3::Eigen)