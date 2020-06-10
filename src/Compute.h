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
	static void Compute5D(const Parameters& params);
	static void Compute6D(const Parameters& params);

private :
	static void ComputeFractal3D(const Parameters& params, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal4Dx0(const Parameters& params, size_t i, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal4Dx1(const Parameters& params, size_t j, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal4Dx2(const Parameters& params, size_t k, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal4Dx3(const Parameters& params, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx0x1(const Parameters& params, size_t i, size_t j, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx0x2(const Parameters& params, size_t i, size_t k, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx0x3(const Parameters& params, size_t i, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx0x4(const Parameters& params, size_t i, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx1x2(const Parameters& params, size_t j, size_t k, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx1x3(const Parameters& params, size_t j, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx1x4(const Parameters& params, size_t j, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx2x3(const Parameters& params, size_t k, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx2x4(const Parameters& params, size_t k, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol);
	static void ComputeFractal5Dx3x4(const Parameters& params, size_t m, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol);
	static unsigned char ComputeMandelBulb4D(float x0, float x1, float x2, float x3, const Parameters& params);
	static unsigned char ComputeMandelBulb5D(float x0, float x1, float x2, float x3, float x4, const Parameters& params);
	static void MeshVolumeAndRender(const Parameters& params, open3d::integration::UniformTSDFVolume& tsdfVol);
}; 


