#pragma once

#include "common/wurst.h"
#include "common/sphericalrepresentation.h"
#include "common/util/tensor.h"

#include "sphericalrepresentation/sphericalfibonacci.h"

template <typename Type>
struct SgLobe
{
	Type eval(const Vec3 & w) const
	{
		assert(Math::IsNormalized(w));
		return mAmplitude * Math::Exp((1.0 - Math::Dot(w, mAxis)) * mSharpness);
	}

	Vec3 mAxis;
	double mSharpness;
	Type mAmplitude;
};

template <typename Type>
struct SgRep : public SphRep<Type>
{
	static void Optimize(SphFiboRep<Type> & map, int numLobes)
	{
		#ifdef WURST_LIBTORCH
			// default pytorch options
			const torch::DeviceType torchDeviceType = torch::kCUDA;
			const torch::ScalarType torchScalarType = torch::kFloat32;
			const torch::TensorOptions option = torch::device(torchDeviceType).dtype(torchScalarType);
			const torch::TensorOptions optionWithGrad = option.requires_grad(true);

			// check map resolution
			int mapResolution = static_cast<int>(map.mData.size());

			// target envmap
			torch::Tensor target = torch::zeros(mapResolution, option);
			for (int i = 0; i < mapResolution; i++)
			{
				target[i] = map.mData[i];
			}

			Fimage<double> image = Util::Tensor::VisualizeTensorFibo(target, Ivec2(1024, 512));
			FimageIo::Save(image, "v.pfm");

			// populate directions for fibonacci sphere
			torch::Tensor fiboX = torch::zeros(mapResolution, option);
			torch::Tensor fiboY = torch::zeros(mapResolution, option);
			torch::Tensor fiboZ = torch::zeros(mapResolution, option);
			for (int i = 0; i < mapResolution; i++)
			{
				Vec3 fiboDir = Mapping::WorldFromSphFiboIndex(i, mapResolution);
				fiboX[i] = fiboDir[0];
				fiboY[i] = fiboDir[1];
				fiboZ[i] = fiboDir[2];
			}

			torch::Tensor unscaledTheta = torch::rand(numLobes, optionWithGrad);
			torch::Tensor unscaledPhi = torch::rand(numLobes, optionWithGrad);

			torch::Tensor unscaledSharpnessTensor = torch::rand(numLobes, optionWithGrad);
			torch::Tensor amplitudeTensor = torch::rand(numLobes, optionWithGrad);

			std::vector<torch::optim::Adam> optimizers;

			// optimizer1 for optimizing angle
			std::vector<torch::Tensor> arguments1;
			arguments1.push_back(unscaledTheta);
			arguments1.push_back(unscaledPhi);
			torch::optim::Adam optimizer1(arguments1, torch::optim::AdamOptions(0.3));
			optimizers.push_back(optimizer1);

			// optimizer2 for optimizing sharpness
			std::vector<torch::Tensor> arguments2;
			arguments2.push_back(unscaledSharpnessTensor);
			torch::optim::Adam optimizer2(arguments2, torch::optim::AdamOptions(1e-2));
			optimizers.push_back(optimizer2);

			// optimizer2 for optimizing amplitude
			std::vector<torch::Tensor> arguments3;
			arguments3.push_back(amplitudeTensor);
			torch::optim::Adam optimizer3(arguments3, torch::optim::AdamOptions(0.05));
			optimizers.push_back(optimizer3);

			for (int iOpt = 0; iOpt < 10000; iOpt++)
			{
				// Sphrerical(theta, phi) --maps to--> Vec3(costheta * sinphi, cosphi, sintheta * sinphi)
				torch::Tensor theta = unscaledTheta * 2.0 * Math::Pi;
				torch::Tensor phi = unscaledPhi * Math::Pi;

				torch::Tensor sinTheta = torch::sin(unscaledTheta);
				torch::Tensor cosTheta = torch::cos(unscaledTheta);
				torch::Tensor sinPhi = torch::sin(unscaledPhi);
				torch::Tensor cosPhi = torch::cos(unscaledPhi);

				torch::Tensor lobeAxisX = cosTheta * sinPhi;
				torch::Tensor lobeAxisY = cosPhi;
				torch::Tensor lobeAxisZ = sinTheta * sinPhi;

				torch::Tensor result;
				for (int i = 0; i < numLobes; i++)
				{
					// compute dot product
					torch::Tensor dot = lobeAxisX[i] * fiboX + lobeAxisY[i] * fiboY + lobeAxisZ[i] * fiboZ;

					// eval sg lobe
					torch::Tensor amplitude = torch::abs(amplitudeTensor[i]); // basically this is clamping
					torch::Tensor sharpness = torch::abs(unscaledSharpnessTensor[i] * 10.0); // basically this is clamping
					torch::Tensor lobe = amplitude * torch::exp((1.0 - dot) * sharpness - sharpness);

					// add the result
					if (i == 0)
						result = lobe;
					else
						result += lobe;
				}
				result /= static_cast<float>(numLobes);
				torch::Tensor renderLoss = (torch::pow(target - result, 2.0) / (result + 0.01)).mean();

				std::cout << renderLoss << std::endl;

				Fimage<double> image = Util::Tensor::VisualizeTensorFibo(result, Ivec2(1024, 512));
				FimageIo::Save(image, "v" + std::to_string(iOpt) + ".hdr");

				try
				{
					optimizers[iOpt % optimizers.size()].zero_grad();
					renderLoss.backward();
					optimizers[iOpt % optimizers.size()].step();
				}
				catch (const std::exception & e)
				{
					std::cout << e.what() << std::endl;
				}
			}
		#endif
	}

	SgRep()
	{
	}

	SgRep(const int numLobes)
	{
	}

	Type eval(const Vec3 & dir) const override
	{
		Type result(0.0);
		for (Uint i = 0; i < mLobes.size(); i++)
		{
			result += mLobes[i].eval(dir);
		}
		return result;
	}

	std::vector<SgLobe<Type>> mLobes;
};
