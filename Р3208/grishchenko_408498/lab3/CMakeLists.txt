cmake_minimum_required(VERSION 3.10.0)
project(cm3 VERSION 0.1.0)

file(GLOB SOURCES *.cpp)

add_executable(cm3 ${SOURCES})
target_link_libraries(cm3 PUBLIC matplot)