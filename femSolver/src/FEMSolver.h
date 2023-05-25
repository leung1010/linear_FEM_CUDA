#ifndef _FEMSOLVER_H_
#define _FEMSOLVER_H_
#include "libFEM.h"
#include "preprocessor.h"
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include "cuFEM.h"

void linearFEMcuSolver(Preprocessor &femData);

Eigen::SparseMatrix<double> AssemblyStiffMat(Preprocessor &femData);

void ApplyConstrainedDOFs(Eigen::SparseMatrix<double>& globalK, const std::vector<ConstrainedDOF>& constrainedDOFs);




#endif //_FEMSOLVER_H_
