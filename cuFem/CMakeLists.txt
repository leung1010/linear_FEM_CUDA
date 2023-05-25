cmake_minimum_required (VERSION 3.5)
find_package(CUDA REQUIRED)

MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")

set(CUDA_64_BIT_DEVICE_CODE ON)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=sm_52")

MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")

CUDA_ADD_LIBRARY(
	cuFEM 
	src/cuFEM.h 
	src/cuFEM.cu
	)

target_include_directories(cuFEM PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(cuFEM ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})
#set_target_properties(cuFEM PROPERTIES LINKER_LANGUAGE CUDA)


