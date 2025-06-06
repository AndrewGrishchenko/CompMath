cmake_minimum_required(VERSION 3.10.0)
project(cm6 VERSION 0.1.0)

include(FetchContent)
FetchContent_Declare(sciplot_content
  GIT_REPOSITORY https://github.com/sciplot/sciplot.git
  GIT_TAG master)

FetchContent_GetProperties(sciplot_content)
if(NOT sciplot_content_POPULATED)
  FetchContent_Populate(sciplot_content)
endif()

include_directories(${sciplot_content_SOURCE_DIR})

file(GLOB SOURCES *.cpp)

add_executable(cm6 ${SOURCES})
target_link_libraries(cm6 PUBLIC matplot)