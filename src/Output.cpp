#include "Output.h"
#include "Parameters.h"
#include <Open3D/Open3D.h>

#include <fstream>
#include <sstream>

template <typename T>
void SwapEnd(T& var)
{
	char* varArray = reinterpret_cast<char*>(&var);
	for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
		std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}

void Output::ExportVtk(const std::string& filePath, const Parameters& params, float step, const std::vector<unsigned char>& voxel)
{
	std::ofstream f;
	f.open(filePath, std::ios::out | std::ios::binary);
	f << "# vtk DataFile Version 2.0" << std::endl;
	f << "Mandelbulb voxel" << std::endl;
	f << "BINARY" << std::endl;
	f << "DATASET STRUCTURED_POINTS" << std::endl;
	f << "DIMENSIONS " << params.N << " " << params.N << " " << params.N << std::endl;
	f << "SPACING " << step << " " << step << " " << step << std::endl;
	f << "ORIGIN " << params.min << " " << params.min << " " << params.min << std::endl;
	f << "POINT_DATA " << params.N * params.N * params.N << std::endl;
	f << "SCALARS volume_scalars float 1" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (size_t v = 0; v < params.N * params.N * params.N; ++v)
	{
		float result = (float)voxel[v] / (float)params.maxIter;
		SwapEnd(result);
		f.write((char*)&result, sizeof(float));
	}
	f.close();
}

void Output::ExportPNG(const std::string& filePath, const open3d::geometry::Image& image)
{
	open3d::io::WriteImageToPNG(filePath, image, 90);
}