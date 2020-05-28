#pragma once

struct Parameters;

namespace open3d {
	namespace integration {
		class UniformTSDFVolume;
	}
}

class Compute {
public :
	static void Compute3D(const Parameters& params);
	static void Compute4D(const Parameters& params);

private :
	static void ComputeFractal3D(const Parameters& params, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal4Dx0(const Parameters& params, size_t i, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal4Dx1(const Parameters& params, size_t j, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal4Dx2(const Parameters& params, size_t k, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal4Dx3(const Parameters& params, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol);
	static unsigned char ComputeMandelBulb4D(float x0, float x1, float x2, float x3, const Parameters& params);
	static void MeshVolumeAndRender(const Parameters& params, size_t id, open3d::integration::UniformTSDFVolume& tsdfVol);
};


