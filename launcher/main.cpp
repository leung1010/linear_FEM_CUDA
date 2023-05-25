#include <iostream>
#include <string>
#include <errno.h>
#include <signal.h>
#include <sys/stat.h>

#include "FEMSolver.h"
#include "preprocessor.h"


int main(int argc, char* argv[])
{
    std::cout << "========GPU-based FEM Simulator============" << std::endl;
    // preprocessor
    // material parameter 
    double E = 1000.0, v = 0.3;
    std::string name = "test";
    std::string project_dir = "/home/liangqx/linearFEMGPU/";
    std::string mesh_dir = project_dir + "mesh/";
    std::string result_dir = project_dir + "result/";
    
	std::string cmd("mkdir -p " + result_dir + name + "/");
	int ret = system(cmd.c_str());
	if (ret)
	{
		std::cout << "create dir error: " << ret << std::endl;
		return -1;
	}
	std::cout << "system cmd create dir successly: " << cmd << std::endl;
    Preprocessor femData;
    femData.parseFiles(mesh_dir, name, E, v);
    femData.printInfo();

    // solver
    linearFEMcuSolver(femData);

    // postprocessor

    return 0;
}