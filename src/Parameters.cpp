#include "Parameters.h"
#include "nlohmann/json.hpp"
#include <fstream>
#include <iostream>

using json = nlohmann::json;

void Parameters::import(const std::string& parametersFileName)
{
	std::ifstream inputFile(parametersFileName);
	json j;
	inputFile >> j;

	dirOutput = j["dirOutput"].get<std::string>();
	std::cout << "Output directory " << dirOutput << std::endl;

	imageFileName = j["imageFileName"].get<std::string>();
	std::cout << "Image name " << imageFileName << std::endl;

	width = j["width"].get<unsigned int>();
	height = j["height"].get<unsigned int>();
	std::cout << "Image width x height " << width << " " << height << std::endl;
	renderingQuality = (RenderingQuality)j["RenderingQuality"].get<int>(); // 0 or 1
	dim = j["dim"].get<unsigned char>(); // 3 to 6
	fixedCoords[0] = j["fixedCoords"][0].get<size_t>();
	fixedCoords[1] = j["fixedCoords"][1].get<size_t>();
	fixedCoords[2] = j["fixedCoords"][2].get<size_t>();
	valueCoords[0] = j["valueCoords"][0].get<size_t>();
	valueCoords[1] = j["valueCoords"][1].get<size_t>();
	valueCoords[2] = j["valueCoords"][2].get<size_t>();
	power = j["power"].get<float>();
	maxIter = j["maxIter"].get<unsigned char>();
	bailout = j["bailout"].get<float>();
	min = j["min"].get<float>();
	extent = j["extent"].get<float>();
	N = j["N"].get<size_t>();
}

bool Parameters::isValid() const 
{
	// only dimensions 3, 4, 5 and 6 
	if ((dim < 3) || (dim > 6))
		return  false;
	if (fixedCoords[0] != size_t(-1))
	{
		if (fixedCoords[0] >= dim)
			return false;
	}
	if (fixedCoords[1] != size_t(-1))
	{
		if (fixedCoords[1] >= dim)
			return false;
	}
	if (fixedCoords[2] != size_t(-1))
	{
		if (fixedCoords[2] >= dim)
			return false;
	}
	if ((dim == 4) && (valueCoords[0] >= N))
		return false;
	if ((dim == 5) && (valueCoords[0] >= N) && (valueCoords[1] >= N))
		return false;
	if ((dim == 6) && (valueCoords[0] >= N) && (valueCoords[1] >= N) && (valueCoords[2] >= N))
		return false;

	return true;
}

std::ostream& operator<<(std::ostream& os, const Parameters& param)
{
	os << "Output directory " << param.dirOutput << std::endl;
	os << "Image name " << param.imageFileName << std::endl;
	os << "Image width x height " << param.width << " " << param.height << std::endl;
	if (param.renderingQuality == RenderingQuality::High)
	{
		os << "Rendering quality High" << std::endl;
	}
	else
	{
		os << "Rendering quality Low" << std::endl;
	}
	os << "Dimension " << static_cast<unsigned>(param.dim) << std::endl;
	if(param.fixedCoords[0] != size_t(-1))
		os << "Fixed coordinates " << param.fixedCoords[0] << std::endl;
	if (param.fixedCoords[1] != size_t(-1))
		os << "Fixed coordinates " << param.fixedCoords[1] << std::endl;
	if (param.fixedCoords[2] != size_t(-1))
		os << "Fixed coordinates " << param.fixedCoords[2] << std::endl;
	if (param.dim == 4)
		os << "Fixed value coordinate " << param.valueCoords[0] << std::endl;
	if (param.dim == 5)
		os << "Fixed value coordinates " << param.valueCoords[0] << " " << param.valueCoords[1] << std::endl;
	if (param.dim == 6)
		os << "Fixed value coordinates " << param.valueCoords[0] << " " << param.valueCoords[1] << " " << param.valueCoords[2] << std::endl;
	os << "Power " << param.power << std::endl;
	os << "Maximum number of iterations " << static_cast<unsigned>(param.maxIter) << std::endl;
	os << "Bailout " << param.bailout << std::endl;
	os << "Min " << param.min << std::endl;
	os << "Extent " << param.extent << std::endl;
	os << "Number sample per axis " << param.N << std::endl;

	return os;
}