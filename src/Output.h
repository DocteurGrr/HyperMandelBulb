#pragma once

#include <string>
#include <vector>

struct Parameters;
namespace open3d {
	namespace geometry {
		class Image;
	}
}

class Output
{
public:
	static void ExportVtk(const std::string& filePath, const Parameters& params, float step, const std::vector<unsigned char>& voxel);
	static void ExportPNG(const std::string& filePath, const open3d::geometry::Image& image);
};