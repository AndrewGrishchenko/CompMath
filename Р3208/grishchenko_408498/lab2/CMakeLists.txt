cmake_minimum_required(VERSION 3.10.0)
project(cm2 VERSION 0.1.0)

include(FetchContent)
FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master)
FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

file(GLOB SOURCES *.cpp)

add_executable(cm2 ${SOURCES})
target_link_libraries(cm2 PUBLIC matplot)