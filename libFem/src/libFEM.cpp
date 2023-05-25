#include "libFEM.h"
#include <iostream>


void Element::localStiffMat(double* NodeX, double* NodeY, Eigen::Matrix<double,3,3> materialD )//
{
    // std::cout << "======= 2D triangle element ==========" << std::endl;
    // std::cout << " compute the local stiff matrix " << std::endl;
    NodesIdx[0]--;NodesIdx[1]--;NodesIdx[2]--;
    Eigen::Vector3d x,y;
    x[0] = NodeX[NodesIdx[0]]; x[1] = NodeX[NodesIdx[1]]; x[2] = NodeX[NodesIdx[2]];
    y[0] = NodeY[NodesIdx[0]]; y[1] = NodeY[NodesIdx[1]]; y[2] = NodeY[NodesIdx[2]];
    // std::cout << NodesIdx[0] <<" " << NodesIdx[1] << " " << NodesIdx[2] << std::endl;
    // std::cout << x << std::endl;
    // std::cout << y << std::endl;
    // std::cout << NodeX[NodesIdx[0]-1] << " " << NodeX[NodesIdx[1]-1] << " " << NodeX[NodesIdx[2]-1] << std::endl;
    // std::cout << NodeY[NodesIdx[0]-1] << " " << NodeY[NodesIdx[1]-1] << " " << NodeY[NodesIdx[2]-1] << std::endl;
    
    // shape function other than {1,x,y}
    Eigen::Matrix3d N;
    N(0,0) = 1.0; N(0,1) = x[0]; N(0,2) = y[0];
    N(1,0) = 1.0; N(1,1) = x[1]; N(1,2) = y[1];
    N(2,0) = 1.0; N(2,1) = x[2]; N(2,2) = y[2];
    // std::cout << N.determinant() << std::endl;
    double TriArea = N.determinant()/2.0;

    Eigen::Matrix3d invN;
    invN = N.inverse();

    for (int i = 0; i < 3; i++)
    {
        B(0, 2*i+0) = invN(1, i);
        B(0, 2*i+1) = 0.0;
        B(1, 2*i+0) = 0.0;
        B(1, 2*i+1) = invN(2, i);
        B(2, 2*i+0) = invN(2, i);
        B(2, 2*i+1) = invN(1, i);
    }
    
    localK = B.transpose()*materialD*B;
    localK *= TriArea;
    // std::cout << localK << std::endl;
}



