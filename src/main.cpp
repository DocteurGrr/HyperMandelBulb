
//#define NOMINMAX

#include "Parameters.h"
#include "Compute.h"

#include <iostream>
#include <chrono>

using namespace std::chrono;

int main(int argc, char** argv)
{
	std::cout << "Hyper MandelBulb computation" << std::endl;

	if (argc != 2)
	{
		std::cout << "Usage HyperMandelBulb pathToConfigFile.json" << std::endl;
		return -1;
	}

	std::string jsonFilePath = std::string(argv[1]);
	std::cout << "Reading " << jsonFilePath << std::endl;
	Parameters params;
	params.import(jsonFilePath);
	if (!params.isValid())
	{
		std::cout << "Invalid config file !" << std::endl;
		return -1;
	}
	std::cout << params << std::endl;

	// start computation
	high_resolution_clock::time_point start = high_resolution_clock::now();
	switch (params.dim)
	{
	case 3 :
		Compute::Compute3D(params);
		break;
	case 4:
		Compute::Compute4D(params);
		break;
	case 5 :
		break;
	case 6 : 
		break;
	default : 
		std::cout << "Computations not implemented for dimension " << params.dim << std::endl;
		return -1;
	}
	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Total Fractal computation " << duration_cast<seconds>(end - start).count() << " s" << std::endl;

	return 0;
}