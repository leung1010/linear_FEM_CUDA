#include <stdio.h>
#include <cuda_runtime.h>
#include <cusolverSp.h>
#include <cusparse_v2.h>
#include <iostream>

#include "cuFEM.h"

void QRwithCusolver(int nnz, int n, double* h_csrValA, double* h_b, double* h_x, int* h_csrRowPtrA, int* h_csrColIndA)
{
    // QR direct way with cusolver
    cusolverSpHandle_t handle;
    cusolverStatus_t solverStatus;

    cusparseStatus_t spStatus;
    cusparseMatDescr_t descrA;

    solverStatus = cusolverSpCreate(&handle);
    spStatus = cusparseCreateMatDescr(&descrA);

    double* d_csrValA, *d_b, *d_x;
    int* d_csrRowPtrA, *d_csrColIndA;
    cudaMalloc((void**)&d_csrValA, nnz * sizeof(double));
    cudaMalloc((void**)&d_b, n * sizeof(double));
    cudaMalloc((void**)&d_x, n * sizeof(double));
    cudaMalloc((void**)&d_csrRowPtrA, (n+1) * sizeof(int));
    cudaMalloc((void**)&d_csrColIndA, n * sizeof(int));

    cudaMemcpy(d_csrValA, h_csrValA, nnz*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, n*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csrRowPtrA, h_csrRowPtrA, (n+1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csrColIndA, h_csrColIndA, n * sizeof(int), cudaMemcpyHostToDevice);

    std::cout << "----------solving with cusolver---------------" << std::endl;
    double tol = 1e-16;
    int reorder = 1;
    int sigularity = 0;
    solverStatus = cusolverSpDcsrlsvqr(handle, n, nnz, descrA, d_csrValA, d_csrRowPtrA, d_csrColIndA, d_b, tol, reorder, d_x, &sigularity);

    std::cout << "----------END solving with cusolver---------------" << std::endl;

    cudaMemcpy(h_x, d_x, n*sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_csrRowPtrA);
    cudaFree(d_csrValA);
    cudaFree(d_csrColIndA);
    cudaFree(d_b);
    cudaFree(d_x);
    cusolverSpDestroy(handle);
    
}
