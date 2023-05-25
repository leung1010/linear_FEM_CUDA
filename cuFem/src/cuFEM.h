#ifndef CUDA_LIB_
#define CUDA_LIB_

void QRwithCusolver(int nnz, int n, double* h_csrValA, double* h_b, double* h_x, int* d_csrRowPtrA, int* d_csrColIndA);


#endif // CUDA_LIB_