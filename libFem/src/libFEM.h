#ifndef _LIBFEM_H_
#define _LIBFEM_H_

#include <Eigen/Dense>

struct Element
{
// private:
//     /* data */
// public:
//     Element(/* args */); // class is OK
//     ~Element();
    void localStiffMat(double* NodeX, double* NodeY, Eigen::Matrix<double,3,3> materialD);
    
    int NodesIdx[3];
    Eigen::Matrix<double, 3, 6> B;
    Eigen::Matrix<double, 6, 6> localK;
};

struct ConstrainedDOF
{
    enum enumDOFs
    {
        Ux=1,Uy,Uz
    };
    int constrainedNode;
    enumDOFs enumdofs;
};







#endif //_LIBFEM_H_