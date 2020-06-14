#include "Compute.h"
#include "Parameters.h"
#include "Render.h"
#include "Output.h"

#include <Open3D/Open3D.h>
#include <Open3D/Integration/UniformTSDFVolume.h>
#include <tbb/tbb.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <sstream>

//#define DEBUG_VOLUME
//#define DEBUG_MESHING
//#define VERBOSE

using namespace std::chrono;
const float M_2PI = 6.28318531f;

void Compute::Compute3D(const Parameters& params)
{
	std::cout << "3D Computation" << std::endl;

	const Eigen::Vector3d origin = Eigen::Vector3d((double)params.min, (double)params.min, (double)params.min);
	double sdf_trunc = -1.5 / (double)params.maxIter;
	open3d::integration::UniformTSDFVolume tsdfVol(
		(double)params.extent,
		(int)params.N,
		sdf_trunc,
		open3d::integration::TSDFVolumeColorType::NoColor,
		origin
	);

	ComputeFractal3D(params, tsdfVol);
  MeshVolumeAndRender(params, tsdfVol);
}

void Compute::Compute4D(const Parameters& params)
{
	std::cout << "4D Computation" << std::endl;

	const Eigen::Vector3d origin = Eigen::Vector3d((double)params.min, (double)params.min, (double)params.min);
	double sdf_trunc = -1.5 / (double)params.maxIter;

	open3d::integration::UniformTSDFVolume tsdfVol(
		(double)params.extent,
		(int)params.N,
		sdf_trunc,
		open3d::integration::TSDFVolumeColorType::NoColor,
		origin
	);

	std::cout << "Computing image " << params.valueCoords[0] << "/" << params.N << std::endl;
	switch (params.fixedCoords[0])
	{
	case 0 :
		ComputeFractal4Dx0(params, params.valueCoords[0], tsdfVol);
		break;
	case 1 :
		ComputeFractal4Dx1(params, params.valueCoords[0], tsdfVol);
		break;
	case 2 :
		ComputeFractal4Dx2(params, params.valueCoords[0], tsdfVol);
		break;
	case 3 :
		ComputeFractal4Dx3(params, params.valueCoords[0], tsdfVol);
		break;
	}

	MeshVolumeAndRender(params, tsdfVol);
	std::cout << "************************************** "<< std::endl;
}

void Compute::Compute5D(const Parameters& params)
{
	std::cout << "5D Computation" << std::endl;

	const Eigen::Vector3d origin = Eigen::Vector3d((double)params.min, (double)params.min, (double)params.min);
	double sdf_trunc = -1.5 / (double)params.maxIter;

	open3d::integration::UniformTSDFVolume tsdfVol(
		(double)params.extent,
		(int)params.N,
		sdf_trunc,
		open3d::integration::TSDFVolumeColorType::NoColor,
		origin
	);

	switch (params.fixedCoords[0])
	{
	case 0:
	{
		switch (params.fixedCoords[1])
		{
		case 1:
			ComputeFractal5Dx0x1(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		case 2:
			ComputeFractal5Dx0x2(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		case 3:
			ComputeFractal5Dx0x3(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		case 4:
			ComputeFractal5Dx0x4(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		}
	}
	break;
	case 1:
	{
		switch (params.fixedCoords[1])
		{
		case 2:
			ComputeFractal5Dx1x2(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		case 3:
			ComputeFractal5Dx1x3(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		case 4:
			ComputeFractal5Dx1x4(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		}
	}
	break;
	case 2:
	{
		switch (params.fixedCoords[1])
		{
		case 3:
			ComputeFractal5Dx2x3(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		case 4:
			ComputeFractal5Dx2x4(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		}
	}
	break;
	case 3:
	{
		switch (params.fixedCoords[1])
		{
		case 4:
			ComputeFractal5Dx3x4(params, params.valueCoords[0], params.valueCoords[1], tsdfVol);
			break;
		}
	}
  break;
	}

	MeshVolumeAndRender(params, tsdfVol);
	std::cout << "************************************** " << std::endl;
}

void Compute::Compute6D(const Parameters& params)
{
	std::cout << "6D Computation" << std::endl;

	const Eigen::Vector3d origin = Eigen::Vector3d((double)params.min, (double)params.min, (double)params.min);
	double sdf_trunc = -1.5 / (double)params.maxIter;

	open3d::integration::UniformTSDFVolume tsdfVol(
		(double)params.extent,
		(int)params.N,
		sdf_trunc,
		open3d::integration::TSDFVolumeColorType::NoColor,
		origin
	);

	switch (params.fixedCoords[0])
	{
	case 0:
	{
		switch (params.fixedCoords[1])
		{
		case 1:
			switch (params.fixedCoords[2])
			{
			case 2:
				ComputeFractal5Dx0x1x2(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 3:
				ComputeFractal5Dx0x1x3(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 4:
				ComputeFractal5Dx0x1x4(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 5:
				ComputeFractal5Dx0x1x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			}
			break;
		case 2:
			switch (params.fixedCoords[2])
			{
			case 3:
				ComputeFractal5Dx0x2x3(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 4:
				ComputeFractal5Dx0x2x4(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 5:
				ComputeFractal5Dx0x2x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			}
			break;
		case 3:
			switch (params.fixedCoords[2])
			{
			case 4:
				ComputeFractal5Dx0x3x4(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 5:
				ComputeFractal5Dx0x3x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			}
			break;
		case 4:
			ComputeFractal5Dx0x4x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
			break;
		}
		break;
	}
	break;
	case 1:
	{
		switch (params.fixedCoords[1])
		{
		case 2:
			switch (params.fixedCoords[2])
			{
			case 3:
				ComputeFractal5Dx1x2x3(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 4:
				ComputeFractal5Dx1x2x4(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 5:
				ComputeFractal5Dx1x2x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			}
			break;
		case 3:
			switch (params.fixedCoords[2])
			{
			case 4:
				ComputeFractal5Dx1x3x4(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 5:
				ComputeFractal5Dx1x3x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			}
			break;
		case 4:
			ComputeFractal5Dx1x4x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
			break;
		}
		break;
	}
	break;
	case 2:
	{
		switch (params.fixedCoords[1])
		{
		case 3:
			switch (params.fixedCoords[2])
			{
			case 4:
				ComputeFractal5Dx2x3x4(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			case 5:
				ComputeFractal5Dx2x3x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
				break;
			}
			break;
		case 4:
			ComputeFractal5Dx2x4x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
			break;
		}
		break;
	}
	break;
	case 3:
	{
		ComputeFractal5Dx3x4x5(params, params.valueCoords[0], params.valueCoords[1], params.valueCoords[2], tsdfVol);
	}
	break;
	}

	MeshVolumeAndRender(params, tsdfVol);
	std::cout << "************************************** " << std::endl;
}

void Compute::ComputeFractal3D(const Parameters& params, open3d::integration::UniformTSDFVolume& tsdfVol)
{
#ifdef DEBUG_VOLUME
	std::vector<unsigned char> voxel(params.N * params.N * params.N);
#endif
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, params.N * params.N * params.N, [&](size_t v)
		{
			size_t i = v / (params.N * params.N);
			size_t t = v % (params.N * params.N);
			size_t j = t / params.N;
			size_t k = t % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;

			unsigned char iter = 0;
			float r = sqrtf(x0 * x0 + x1 * x1 + x2 * x2);

			while ((r < params.bailout) && (iter < params.maxIter))
			{
				// convert to polar coordinates
				float theta = acosf(x2 / r);
				float phi = atan2f(x1, x0);

				// scale and rotate the point
				float x2r = powf(r, params.power);
				theta *= params.power;
				phi *= params.power;

				// convert back to cartesian coordinates
				x0 += x2r * sin(theta) * cos(phi);
				x1 += x2r * sin(theta) * sin(phi);
				x2 += x2r * cos(theta);

				// prepare next iteration
				r = sqrtf(x0 * x0 + x1 * x1 + x2 * x2);
				iter++;
			}

#ifdef DEBUG_VOLUME
			voxel[v] = iter;
#endif
			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)k));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;

	// save to vtk file
#ifdef DEBUG_VOLUME
	std::string filepath = params.dirOutput + params.imageFileName + std::string(".vtk");
	Output::ExportVtk(filepath, params, step, voxel);
	std::cout << "Voxel file " << filepath << " saved !" << std::endl;
#endif
}

void Compute::ComputeFractal4Dx0(const Parameters& params, size_t i, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t m = v / N2;
			size_t t = v % N2;
			size_t k = t / params.N;
			size_t j = t % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb4D(x0, x1, x2, x3, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)k, (int)m));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal4Dx1(const Parameters& params, size_t j, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t m = v / N2;
			size_t t = v % N2;
			size_t k = t / params.N;
			size_t i = t % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb4D(x0, x1, x2, x3, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)k, (int)m));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal4Dx2(const Parameters& params, size_t k, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t m = v / N2;
			size_t t = v % N2;
			size_t j = t / params.N;
			size_t i = t % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb4D(x0, x1, x2, x3, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)m));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal4Dx3(const Parameters& params, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t k = v / N2;
			size_t t = v % N2;
			size_t j = t / params.N;
			size_t i = t % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb4D(x0, x1, x2, x3, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)k));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x1(const Parameters& params, size_t i, size_t j, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t m = r / params.N;
			size_t k = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)k, (int)m, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x2(const Parameters& params, size_t i, size_t k, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t m = r / params.N;
			size_t j = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)m, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x3(const Parameters& params, size_t i, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t j = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)k, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x4(const Parameters& params, size_t i, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t m = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t j = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)k, (int)m));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx1x2(const Parameters& params, size_t j, size_t k, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t m = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)m, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx1x3(const Parameters& params, size_t j, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)k, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx1x4(const Parameters& params, size_t j, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t m = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)k, (int)m));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx2x3(const Parameters& params, size_t k, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t j = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx2x4(const Parameters& params, size_t k, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t m = v / N2;
			size_t r = v % N2;
			size_t j = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)m));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx3x4(const Parameters& params, size_t m, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t k = v / N2;
			size_t r = v % N2;
			size_t j = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb5D(x0, x1, x2, x3, x4, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)k));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x1x2(const Parameters& params, size_t i, size_t j, size_t k, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t p = r / params.N;
			size_t m = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)m, (int)p, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x1x3(const Parameters& params, size_t i, size_t j, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t p = r / params.N;
			size_t k = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)k, (int)p, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x1x4(const Parameters& params, size_t i, size_t j, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t m = r / params.N;
			size_t k = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)k, (int)m, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x1x5(const Parameters& params, size_t i, size_t j, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t m = r / params.N;
			size_t k = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)k, (int)m, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x2x3(const Parameters& params, size_t i, size_t k, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t p = r / params.N;
			size_t j = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)p, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x2x4(const Parameters& params, size_t i, size_t k, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t m = r / params.N;
			size_t j = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)m, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x2x5(const Parameters& params, size_t i, size_t k, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t m = r / params.N;
			size_t j = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)m, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x3x4(const Parameters& params, size_t i, size_t m, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t j = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)k, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x3x5(const Parameters& params, size_t i, size_t m, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t j = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)k, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx0x4x5(const Parameters& params, size_t i, size_t p, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t m = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t j = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)j, (int)k, (int)m));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx1x2x3(const Parameters& params, size_t j, size_t k, size_t m, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t p = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)p, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx1x2x4(const Parameters& params, size_t j, size_t k, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t m = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)m, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx1x2x5(const Parameters& params, size_t j, size_t k, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t m = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)m, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx1x3x4(const Parameters& params, size_t j, size_t m, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)k, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx1x3x5(const Parameters& params, size_t j, size_t m, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)k, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx1x4x5(const Parameters& params, size_t j, size_t p, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t m = v / N2;
			size_t r = v % N2;
			size_t k = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)k, (int)m));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx2x3x4(const Parameters& params, size_t k, size_t m, size_t p, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t t = v / N2;
			size_t r = v % N2;
			size_t j = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)t));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx2x3x5(const Parameters& params, size_t k, size_t m, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t p = v / N2;
			size_t r = v % N2;
			size_t j = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)p));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx2x4x5(const Parameters& params, size_t k, size_t p, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t m = v / N2;
			size_t r = v % N2;
			size_t j = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)m));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}

void Compute::ComputeFractal5Dx3x4x5(const Parameters& params, size_t m, size_t p, size_t t, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	size_t N3 = params.N * params.N * params.N;
	size_t N2 = params.N * params.N;
	float step = params.extent / (float)params.N;

	high_resolution_clock::time_point start = high_resolution_clock::now();

	tbb::parallel_for((size_t)0, N3, [&](size_t v)
		{
			size_t k = v / N2;
			size_t r = v % N2;
			size_t j = r / params.N;
			size_t i = r % params.N;

			float x0 = params.min + ((float)i + 0.5f) * step;
			float x1 = params.min + ((float)j + 0.5f) * step;
			float x2 = params.min + ((float)k + 0.5f) * step;
			float x3 = params.min + ((float)m + 0.5f) * step;
			float x4 = params.min + ((float)p + 0.5f) * step;
			float x5 = params.min + ((float)t + 0.5f) * step;

			unsigned char iter = ComputeMandelBulb6D(x0, x1, x2, x3, x4, x5, params);

			open3d::geometry::TSDFVoxel voxel(Eigen::Vector3i((int)i, (int)j, (int)k));
			voxel.tsdf_ = ((float)iter / (float)params.maxIter) - 1.f;
			voxel.weight_ = 1.f;
			tsdfVol.voxels_[v] = voxel;
		});

	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "Fractal computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}


unsigned char Compute::ComputeMandelBulb4D(float x0, float x1, float x2, float x3, const Parameters &params)
{
	unsigned char iter = 0;

	float r2 = x2 * x2 + x3 * x3;
	float r3 = x1 * x1 + r2;
	float r4 = x0 * x0 + r3;
	float r = sqrtf(r4);

	while ((r < params.bailout) && (iter < params.maxIter))
	{
		// convert to polar coordinates
		// See wikipedia articles on n-sphere
		float phi1 = acosf(x0 / r);
		float phi2 = acosf(x1 / sqrtf(r3));
		float phi3 = (x3 >= 0.f) ? acosf(x2 / sqrtf(r2)) : M_2PI - acosf(x2 / sqrtf(r2));

		// scale and rotate the point
		float zr = powf(r, params.power);
		phi1 *= params.power;
		phi2 *= params.power;
		phi3 *= params.power;

		// convert back to cartesian coordinates
		float cosphi1 = cosf(phi1);
		float sinphi1 = sinf(phi1);
		float cosphi2 = cosf(phi2);
		float sinphi2 = sinf(phi2);
		float cosphi3 = cosf(phi3);
		float sinphi3 = sinf(phi3);

		x0 += zr * cosphi1;
		x1 += zr * sinphi1 * cosphi2;
		x2 += zr * sinphi1 * sinphi2 * cosphi3;
		x3 += zr * sinphi1 * sinphi2 * sinphi3;

		// prepare next iteration
		r2 = x2 * x2 + x3 * x3;
		r3 = x1 * x1 + r2;
		r4 = x0 * x0 + r3;
		r = sqrtf(r4);
		iter++;
	}

	return iter;
}

unsigned char Compute::ComputeMandelBulb5D(float x0, float x1, float x2, float x3, float x4, const Parameters& params)
{
	unsigned char iter = 0;

	float r1 = x3 * x3 + x4 * x4;
	float r2 = x2 * x2 + r1;
	float r3 = x1 * x1 + r2;
	float r4 = x0 * x0 + r3;
	float r = sqrtf(r4);

	while ((r < params.bailout) && (iter < params.maxIter))
	{
		// convert to polar coordinates
		// See wikipedia articles on n-sphere
		float phi1 = acosf(x0 / r);
		float phi2 = acosf(x1 / sqrtf(r3));
		float phi3 = acosf(x2 / sqrtf(r2));
		float phi4 = (x4 >= 0.f) ? acosf(x3 / sqrtf(r1)) : M_2PI - acosf(x3 / sqrtf(r1));

		// scale and rotate the point
		float zr = powf(r, params.power);
		phi1 *= params.power;
		phi2 *= params.power;
		phi3 *= params.power;
		phi4 *= params.power;

		// convert back to cartesian coordinates
		float cosphi1 = cosf(phi1);
		float sinphi1 = sinf(phi1);
		float cosphi2 = cosf(phi2);
		float sinphi2 = sinf(phi2);
		float cosphi3 = cosf(phi3);
		float sinphi3 = sinf(phi3);
		float cosphi4 = cosf(phi4);
		float sinphi4 = sinf(phi4);

		x0 += zr * cosphi1;
		x1 += zr * sinphi1 * cosphi2;
		x2 += zr * sinphi1 * sinphi2 * cosphi3;
		x3 += zr * sinphi1 * sinphi2 * sinphi3 * cosphi4;
		x4 += zr * sinphi1 * sinphi2 * sinphi3 * sinphi4;

		// prepare next iteration
		r1 = x3 * x3 + x4 * x4;
	  r2 = x2 * x2 + r1;
	  r3 = x1 * x1 + r2;
	  r4 = x0 * x0 + r3;
	  r = sqrtf(r4);
		iter++;
	}

	return iter;
}

unsigned char Compute::ComputeMandelBulb6D(float x0, float x1, float x2, float x3, float x4, float x5, const Parameters& params)
{
	unsigned char iter = 0;

	float r0 = x4 * x4 + x5 * x5;
	float r1 = x3 * x3 + r0;
	float r2 = x2 * x2 + r1;
	float r3 = x1 * x1 + r2;
	float r4 = x0 * x0 + r3;
	float r = sqrtf(r4);

	while ((r < params.bailout) && (iter < params.maxIter))
	{
		// convert to polar coordinates
		// See wikipedia articles on n-sphere
		float phi1 = acosf(x0 / r);
		float phi2 = acosf(x1 / sqrtf(r3));
		float phi3 = acosf(x2 / sqrtf(r2));
		float phi4 = acosf(x3 / sqrtf(r1));
		float phi5 = (x5 >= 0.f) ? acosf(x4 / sqrtf(r0)) : M_2PI - acosf(x4 / sqrtf(r0));

		// scale and rotate the point
		float zr = powf(r, params.power);
		phi1 *= params.power;
		phi2 *= params.power;
		phi3 *= params.power;
		phi4 *= params.power;
		phi5 *= params.power;

		// convert back to cartesian coordinates
		float cosphi1 = cosf(phi1);
		float sinphi1 = sinf(phi1);
		float cosphi2 = cosf(phi2);
		float sinphi2 = sinf(phi2);
		float cosphi3 = cosf(phi3);
		float sinphi3 = sinf(phi3);
		float cosphi4 = cosf(phi4);
		float sinphi4 = sinf(phi4);
		float cosphi5 = cosf(phi5);
		float sinphi5 = sinf(phi5);

		x0 += zr * cosphi1;
		x1 += zr * sinphi1 * cosphi2;
		x2 += zr * sinphi1 * sinphi2 * cosphi3;
		x3 += zr * sinphi1 * sinphi2 * sinphi3 * cosphi4;
		x4 += zr * sinphi1 * sinphi2 * sinphi3 * sinphi4 * cosphi5;
		x5 += zr * sinphi1 * sinphi2 * sinphi3 * sinphi4 * sinphi5;

		// prepare next iteration
		r0 = x4 * x4 + x5 * x5;
		r1 = x3 * x3 + r0;
		r2 = x2 * x2 + r1;
		r3 = x1 * x1 + r2;
		r4 = x0 * x0 + r3;
		r = sqrtf(r4);
		iter++;
	}

	return iter;
}


void Compute::MeshVolumeAndRender(const Parameters& params, open3d::integration::UniformTSDFVolume& tsdfVol)
{
	high_resolution_clock::time_point start = high_resolution_clock::now();
	std::shared_ptr<open3d::geometry::TriangleMesh> mesh = tsdfVol.ExtractTriangleMesh();
	mesh->ComputeVertexNormals();
	high_resolution_clock::time_point end = high_resolution_clock::now();
	std::cout << "tsdf marching cube computation " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
	std::cout << "Mesh has " << mesh->triangles_.size() << " triangles and " << mesh->vertices_.size() << " vertex" << std::endl;

#ifdef VERBOSE
	std::cout << "Mesh bounds " << std::endl;
	Eigen::Vector3d minBound = mesh->GetMinBound();
	Eigen::Vector3d maxBound = mesh->GetMaxBound();
	std::cout << minBound[0] << " " << minBound[1] << " " << minBound[2] << std::endl;
	std::cout << maxBound[0] << " " << maxBound[1] << " " << maxBound[2] << std::endl;
#endif

	tsdfVol.Reset();

#ifdef DEBUG_MESHING
	// save	mesh to ply file
	std::stringstream ssp;
	ssp << params.dirOutput << params.imageFileName << "_" << id <<  std::string(".ply");
	std::string meshFilePath = ssp.str();
	open3d::io::WriteTriangleMeshToPLY(meshFilePath, *mesh, false, true, true, true, false, true);
	std::cout << "Mesh file saved to " << meshFilePath << std::endl;
#endif

	open3d::geometry::Image image;
	image.Prepare(params.width, params.height, 3, 1);
	start = high_resolution_clock::now();
	Render::render(mesh, params, image);
	end = high_resolution_clock::now();
	std::cout << "Rendering image in " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
	std::stringstream ss;
	ss << params.dirOutput << params.imageFileName << std::string(".png");
	std::string imageFilePath = ss.str();
	Output::ExportPNG(imageFilePath, image);
	std::cout << "Saved image " << imageFilePath << std::endl;
}