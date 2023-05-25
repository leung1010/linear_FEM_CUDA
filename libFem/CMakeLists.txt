cmake_minimum_required (VERSION 3.10)



add_library(
	libFem 
	src/libFEM.h 
	src/libFEM.cpp
	)

target_include_directories(libFem PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(libFem)