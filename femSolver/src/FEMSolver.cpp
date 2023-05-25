#include "FEMSolver.h"
#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/IterativeLinearSolvers>

//#include "directQRWithCUDA.cuh"
//#include "libFEM.h"

void linearFEMcuSolver(Preprocessor &femData)
{
    // compute the local stiffness matrix
    std::cout << "--------------------- compute the local stiffness matrix --------------------" << std::endl;
    for (std::vector<Element>::iterator iele = femData.elements.begin() ; iele != femData.elements.end(); iele++)
    {
        iele->localStiffMat(femData.coordX, femData.coordY, femData.materialD);
    }
    
    // assembly the global stiffness matrix
    std::cout << "--------------------- assembly the global stiffness matrix --------------------" << std::endl;
    Eigen::SparseMatrix<double> GlobalStiffMat =  AssemblyStiffMat(femData);
    // std::cout << GlobalStiffMat << std::endl;
    GlobalStiffMat.makeCompressed();
    // GlobalStiffMat.isCompressed();
    //std::cout << GlobalStiffMat.isCompressed() << std::endl;
    // std::cout << "sparse for cusolver:    " << GlobalStiffMat.innerIndexPtr() << std::endl;
    // std::cout << "sparse for cusolver:    " << GlobalStiffMat.outerIndexPtr() << std::endl;
    // std::cout << "sparse for cusolver:    " << GlobalStiffMat.valuePtr() << std::endl;
    // std::cout << "sparse for cusolver:    " << GlobalStiffMat.nonZeros() << std::endl;
    // std::cout << "sparse for cusolver:    " << GlobalStiffMat.rows() << std::endl;

    // apply the boundary contion(using constrained dofs)
    std::cout << "--------------------- apply the boundary contion (using constrained dofs) --------------------" << std::endl;
    ApplyConstrainedDOFs(GlobalStiffMat, femData.constrainedDOFs);
    
    std::cout << "---------------------------------- apply the loads --------------------------------" << std::endl;
    Eigen::VectorXd rhs(32);
    double curhs[32]{0.0};
    for (int i = 0; i < 32; i++)
    {
        rhs[i] = femData.loads[i];
        curhs[i] = femData.loads[i];
    }
    //std::cout << rhs << std::endl;


    // solve the linear system A*x = b
    std::cout << "--------------------- solve the linear system A*x = b --------------------" << std::endl;
    Eigen::VectorXd disp(32);
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
    solver.compute(GlobalStiffMat);
    if (solver.info()!=Eigen::Success)
    {
        std::cout << "Eigen ConjugateGradient failed!" << std::endl;
        return;
    }
    disp = solver.solve(rhs);
    std::cout << "Eigen solved displacement is: \n"<< disp << std::endl;

    // solver the linear system with cusolver
    double cudisp[32]{0.0};
    QRwithCusolver(GlobalStiffMat.nonZeros(), GlobalStiffMat.cols(), GlobalStiffMat.valuePtr(), 
                    curhs, cudisp, GlobalStiffMat.outerIndexPtr(), GlobalStiffMat.innerIndexPtr());

    std::cout << "GPU solved displacement is: \n"<< disp << std::endl;
}

Eigen::SparseMatrix<double> AssemblyStiffMat(Preprocessor &femData)
{
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::SparseMatrix<double> SparseGlobalK(32, 32);
    for (std::vector<Element>::iterator iele = femData.elements.begin(); iele != femData.elements.end(); iele++)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int ilocalk = 0; ilocalk < 2; ilocalk++)
                {
                    for (int jlocalk = 0; jlocalk < 2; jlocalk++)
                    {
                        // read the value from the Dense local Stiffness matrix
                        double value = iele->localK(2*i+ilocalk, 2*j+jlocalk);
                        
                        if (value != 0.0)
                        {
                            triplets.push_back(Eigen::Triplet<double>(2*(iele->NodesIdx[i])+ilocalk, 2*(iele->NodesIdx[j])+jlocalk, value));
                        }
                    }
                }
            }
        }
    }
    SparseGlobalK.setFromTriplets(triplets.begin(), triplets.end());

    return SparseGlobalK;
}

void ApplyConstrainedDOFs(Eigen::SparseMatrix<double>& globalK , const std::vector<ConstrainedDOF>& constrainedDOFs)
{
    std::cout << "------------ apply the constraints ----------------" << std::endl;
    std::vector<int> constrainedDofIds; // store the indices of the constrainedDOFs 

    for (std::vector<ConstrainedDOF>::const_iterator icons = constrainedDOFs.begin(); icons != constrainedDOFs.end(); icons++)
    {
        if (icons->enumdofs & ConstrainedDOF::Ux)
        {
            constrainedDofIds.push_back(2*(icons->constrainedNode-1)+0); // -1 is for the constrainedNode 1-indexed
        }
        if (icons->enumdofs & ConstrainedDOF::Uy)
        {
            constrainedDofIds.push_back(2*(icons->constrainedNode-1)+1);
        }        
    }

    int tmp=0;
    for (int i = 0; i < constrainedDofIds.size(); i++)
    {
        tmp = constrainedDofIds[i];// avoid the core dumped
        for (int j = 0; j < 32; j++)
        {
            globalK.coeffRef(tmp,j) = 0.0;
            globalK.coeffRef(j,tmp) = 0.0;
        }
        globalK.coeffRef(tmp,tmp) = 1.0;
    }
}
    