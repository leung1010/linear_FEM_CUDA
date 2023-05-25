#aux_source_directory(. DIR_LIB_SRCS)

#add_library(femSolver STATIC ${DIR_LIB_SRCS})

cmake_minimum_required (VERSION 3.10)


add_library(
	FEMSolver 
	src/FEMSolver.h 
	src/FEMSolver.cpp
	)

target_include_directories(FEMSolver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(FEMSolver preprocessor libFem cuFEM)



