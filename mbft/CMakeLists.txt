cmake_minimum_required(VERSION 3.14)

project(mbft)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src source_list)

add_executable(mbft ${source_list})