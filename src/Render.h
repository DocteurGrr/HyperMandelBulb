#pragma once

#include <memory>

namespace open3d {
	namespace geometry {
		class TriangleMesh;
		class Image;
	}
}
struct Parameters;

class Render {
public :
	static void render(std::shared_ptr<open3d::geometry::TriangleMesh>& mesh, const Parameters& params, open3d::geometry::Image& imageOpen3D);

};

