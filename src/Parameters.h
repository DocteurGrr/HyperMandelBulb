#pragma once

#include <string>

enum class RenderingQuality : int {
	Low = 0,
	High = 1
};

struct Parameters {
	void import(const std::string& parametersFileName);
	bool isValid() const;
	friend std::ostream& operator<<(std::ostream& os, const Parameters& param);

	// Image directory path
	std::string dirOutput;
	// Images root filename + _i .png
	std::string imageFileRootName;
	// Image width and height
	unsigned int width;
	unsigned int height;
	// See enum dfined before Low = 0, High = 1
	RenderingQuality renderingQuality;
	// Fractal number of parameters = between 3 and 6 (3D/4D/5D/6D)
	unsigned char dim; 
	// Which coordinates are fixed ? In order to create 3D volume
	// Fractal 3D = none
	// Fractal 4D = 0 <= fixedCoords[0] <= 3 
	// Fractal 5D = 0 <= fixedCoords[0] <= 4, 0 <= fixedCoords[1] <= 4 
	// Fractal 6D = 0 <= fixedCoords[0] <= 5, 0 <= fixedCoords[1] <= 5, 0 <= fixedCoords[2] <= 5 
	size_t fixedCoords[3];
	// Mandelbuld power (8 for the most famous one)
	float power;
	// Max number of iterations - 1rst criterion for iteration stop
	unsigned char maxIter;
	// 2nd criterion for iteration stop
	float bailout;
	// 3D box min along each axis
	float min;
	// 3D box extent alons each axis 
	float extent;
	// Number of cells along one axis (typical value = 256, 512, 1024). 3D volume has N^3 voxels
	size_t N;
};

std::ostream& operator<<(std::ostream& os, const Parameters& param);