#include "preprocessor.h"
#include <iostream>
#include <fstream>

Preprocessor::Preprocessor(/* args */)
{
}

Preprocessor::~Preprocessor()
{
}

void Preprocessor::parseFiles(std::string dir, std::string name, double E, double v)
{
    std::cout << "------ parsing files... -------------" << std::endl;

    std::fstream nodesFile, elementsFile, loadsFile, constrainedDOFsFile;

    nodesFile.open(dir+name+"/nodes.txt", std::fstream::in);
    elementsFile.open(dir+name+"/elements.txt", std::fstream::in);
    constrainedDOFsFile.open(dir+name+"/constrainedDOFs.txt", std::fstream::in);
    loadsFile.open(dir+name+"/loads.txt", std::fstream::in);

    // read the number of the nodes, elements
    nodesFile >> numNodes;
    elementsFile >> numElements;
    constrainedDOFsFile >> numConstrainedDOFs;
    loadsFile >> numLoads;
    
    // plane stress material matrix
    materialD << 1.0, v, 0.0,
                 v, 1.0, 0.0,
                 0.0, 0.0, (1-v)/2;
    materialD *= E/(1.0-v*v);

    // read the coodinates from the nodes.txt
    // std::cout << numNodes << std::endl;
    for (int i = 0; i < numNodes; i++)
    {
        nodesFile >> coordX[i] >> coordY[i];
        // std::cout << coordX[i] << " " << coordY[i] << std::endl;
    }

    // std::cout << numElements << std::endl;
    for (int i = 0; i < numElements; i++)
    {
        Element element;
        elementsFile >> element.NodesIdx[0] >> element.NodesIdx[1] >> element.NodesIdx[2];
        elements.push_back(element);
        // std::cout << element.NodesIdx[0] << " " << element.NodesIdx[1] << " " << element.NodesIdx[2] << std::endl;
    }

    for (int i = 0; i < numConstrainedDOFs; i++)
    {
        ConstrainedDOF constrainedDOF;
        int constrainedType; // for read the txt file
        constrainedDOFsFile >> constrainedDOF.constrainedNode >> constrainedType;
        constrainedDOF.enumdofs = static_cast<ConstrainedDOF::enumDOFs>(constrainedType);
        constrainedDOFs.push_back(constrainedDOF);
        // std::cout << constrainedDOF.constrainedNode << " " << constrainedType << std::endl;
    }

    for (int i = 0; i < numLoads; i++)
    {
        int loadNode;
        double loadX, loadY;
        loadsFile >> loadNode >> loadX >> loadY;
        loadNode--; //for the load 1-indexed
        loads[2*loadNode+0] = loadX;
        loads[2*loadNode+1] = loadY; 
        // std::cout << loads[2*loadNode+0] << loads[2*loadNode+1] << std::endl;
    }
}

void Preprocessor::printInfo()
{
    std::cout << "================== Data has been prepared for the FEM ============" << std::endl;
}
