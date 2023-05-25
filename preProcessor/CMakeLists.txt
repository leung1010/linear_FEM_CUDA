cmake_minimum_required (VERSION 3.10)


add_library(
	preprocessor 
	src/preprocessor.h 
	src/preprocessor.cpp
	)

target_include_directories(preprocessor PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(preprocessor libFem)