#ifndef _PREPROCESSOR_
#define _PREPROCESSOR_

#include <string>
#include <Eigen/Dense>
#include <vector>
#include "libFEM.h"



class Preprocessor
{
private:
    /* data */
public:
    Preprocessor(/* args */);
    ~Preprocessor();

    void parseFiles(std::string dir, std::string name, double E, double v);
    void printInfo();

    int numNodes;
    int numElements;
    int numConstrainedDOFs;
    int numLoads;
    Eigen::Matrix3d materialD;

    double coordX[16]; // how to store the coodinates, data structure
    double coordY[16];
    std::vector<Element> elements;
    std::vector<ConstrainedDOF> constrainedDOFs;
    double loads[32]{0.0};
};




#endif // _PREPROCESSOR_