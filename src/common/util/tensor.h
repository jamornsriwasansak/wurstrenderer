#pragma once

#include "common/wurst.h"
#include "common/floatimage.h"

#ifdef WURST_LIBTORCH
namespace Util
{
	struct Tensor
	{
		static void SaveFloat32Tensor2d(const filesystem::path & filename, torch::Tensor & t)
		{
			std::ofstream os(filename, std::ios::binary | std::ios::out);
			if (!os.is_open()) throw std::runtime_error("can't save file since it can't be opened");

			torch::Tensor u = t.to(torch::kCPU, torch::kFloat32);

			auto sizes = t.sizes();
			int64_t dim = (int64_t)(sizes.size());
			int64_t numData = 1;
			for (Uint iDim = 0; iDim < sizes.size(); iDim++) numData *= sizes[iDim];

			// write dim
			os.write(reinterpret_cast<char *>(&dim), sizeof(int64_t));

			// write size
			os.write(reinterpret_cast<const char *>(&sizes[0]), sizeof(int64_t) * dim);

			// write data
			os.write(reinterpret_cast<char *>(u.data_ptr()), sizeof(float) * numData);

			// finish writing
			os.close();
		}

		static torch::Tensor LoadFloat32Tensor2d(const filesystem::path & filename, torch::TensorOptions option)
		{
			std::ifstream is(filename, std::ios::binary | std::ios::out);
			if (!is.is_open()) throw std::runtime_error("can't load file since it can't be opened");

			// read dim
			int64_t dim = 0;
			{
				size_t extracted = is.read(reinterpret_cast<char*>(&dim), sizeof(int64_t)).gcount();
				if (extracted != sizeof(int64_t)) throw std::runtime_error("extracted size mismatch");
			}

			// read sizes
			std::vector<int64_t> sizes(dim);
			{
				size_t extracted = is.read(reinterpret_cast<char*>(&sizes[0]), sizeof(int64_t) * dim).gcount();
				if (extracted != sizeof(int64_t) * dim) throw std::runtime_error("extracted size mismatch");
			}

			int64_t numData = 1;
			torch::Tensor result = torch::zeros(sizes).to(torch::kFloat32);
			for (Uint iDim = 0; iDim < sizes.size(); iDim++) numData *= sizes[iDim];
			{
				size_t extracted = is.read(reinterpret_cast<char*>(result.data_ptr()), sizeof(float) * numData).gcount();
				if (extracted != sizeof(float) * numData) throw std::runtime_error("extracted size mismatch");
			}

			return result.to(option);
		}

		static std::pair<torch::Tensor, torch::Tensor> Pca(const torch::Tensor & data2)
		{
			torch::Tensor data = data2.to(torch::kCPU);

			// extract pca from
			torch::Tensor mean = torch::mean(data, 0);
			torch::Tensor centered = data - mean.expand_as(data);

			// find svd and return v
			std::tuple<torch::Tensor, torch::Tensor, torch::Tensor> usv = torch::svd(centered);
			torch::Tensor v = std::get<2>(usv);

			return make_pair(mean.to(data2.device()), torch::t(v.to(data2.device())));
		}

		static torch::Tensor LeftPseudoInverse(torch::Tensor a)
		{
			return torch::mm(torch::mm(a.t(), a).inverse(), a.t());
		}

		static torch::Tensor RightPseudoInverse(torch::Tensor a)
		{
			return torch::mm(a.t(), torch::mm(a, a.t()).inverse());
		}

		template <typename Scalar>
		__forceinline static Scalar ValueFromTensor(const torch::Tensor & t)
		{
			if (t.scalar_type() == torch::kFloat32)
			{
				return static_cast<Scalar>(*(float*)(t.to(torch::kCPU).data_ptr()));
			}
			assert(false); // unhandle
			return Scalar(0);
		}

		template <typename Scalar>
		__forceinline static std::vector<Scalar> Values1dFromTensor(const torch::Tensor & t, const Uint dim = 0)
		{
			Uint size = t.size(dim);
			std::vector<Scalar> result(size);
			if (t.scalar_type() == torch::kFloat32)
			{
				torch::Tensor tcpu = t.to(torch::kCPU);
				float * tcpuPtr = (float*)(t.to(torch::kCPU).data_ptr());

				for (Uint i = 0; i < size; i++) { result[i] = static_cast<Scalar>(tcpuPtr[i]); }
				return result;
			}
			assert(false);
		}

		static Fimage<double> VisualizeTensorFibo(const torch::Tensor & t, const Ivec2 & size)
		{
			torch::Tensor temp = t.to(torch::kCPU);

			Fimage<double> result(size);
			const Uvec2 usize(size);
			const Vec2 dsize(size);
			const int basisResolution = static_cast<int>(t.size(0));
			for (int y = 0; y < usize[1]; y++)
				for (int x = 0; x < usize[0]; x++)
				{
					Vec2 xy(static_cast<double>(x), static_cast<double>(y));
					Vec2 uv = xy / dsize;
					Vec3 direction = Mapping::WorldFromPanorama(uv);
					int fiboIndex = Mapping::SphFiboIndexFromWorld(direction, basisResolution);
					result.at(x, y) = Util::Tensor::ValueFromTensor<double>(temp[fiboIndex]);
				}
			return result;
		}

		static Fimage<double> VisualizeTensor2d(const torch::Tensor & t2)
		{
			const torch::Tensor t = t2.to(torch::kCPU);
			Fimage<double> result(Ivec2(static_cast<int>(t.size(0)), static_cast<int>(t.size(1))));
			for (int r = 0; r < t.size(0); r++)
				for (int c = 0; c < t.size(1); c++)
				{
					result.at(c, r) = (double)ValueFromTensor<float>(t[r][c]);
				}
			return result;
		}
	};
}
#endif
