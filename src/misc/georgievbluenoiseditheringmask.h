#pragma once

#include "common/wurst.h"
#include "common/floatimage.h"
#include "sampler/randomsampler.h"

// my stupidly bruteforce implementation
struct GeorgievBlueNoiseDitheringMask
{
	static double ComputeEnergy(const std::vector<double> & mask, const Ivec2 & imageSize, const int dimension, const double sigmaI = 2.1, const double sigmaS = 1.0)
	{
		double energy = 0.0;
		assert(mask.size() == imageSize[0] * imageSize[1] * dimension);
		Ibound2 bound(Ivec2(0), imageSize);
		const double dimBy2 = double(dimension) * 0.5;
		for (const Ivec2 & i : bound)
			for (const Ivec2 & j : bound)
			{
				if (i != j)
				{
					// compute left term - Distance between pixels
					const int pixelDistance = Math::Distance2(i, j);
					const double l = double(pixelDistance) / (sigmaI * sigmaI);

					const int pOffset = Math::Index(bound, i) * dimension;
					const int qOffset = Math::Index(bound, j) * dimension;

					// compute right term -	Difference between two samples
					double sum = 0.0;
					for (int iDim = 0; iDim < dimension; iDim++)
					{
						const double a = mask[pOffset + iDim];
						const double b = mask[qOffset + iDim];
						sum += Math::Square(a - b);
					}
					const double r = std::pow(std::sqrt(sum), dimBy2) / (sigmaS * sigmaS);
					energy += std::exp(-l -r);
				}
			}
		return energy;
	}

	static void SwapPixel(std::vector<double> * mask, const Ivec2 & imageSize, const int dimension, const Ivec2 & p1, const Ivec2 & p2)
	{
		Ibound2 bound(Ivec2(0), imageSize);
		const int p1Offset = Math::Index(bound, p1) * dimension;
		const int p2Offset = Math::Index(bound, p2) * dimension;

		for (int iDim = 0; iDim < dimension; iDim++)
		{
			std::swap(mask->at(iDim + p1Offset), mask->at(iDim + p2Offset));
		}
	}

	static void SimulateAnnealing(const Ivec2 & imageSize, const int dimension = 1)
	{
		std::vector<double> currentMask(imageSize[0] * imageSize[1] * dimension, 0.0);

		// initialize with white noise with my favourite number
		RandomSampler sampler(1212312121);
		Ibound2 bound(Ivec2(0), imageSize);
		for (const Ivec2 & p : bound)
		{
			int index = Math::Index(bound, p);
			for (int iDim = 0; iDim < dimension; iDim++) { currentMask[index] = sampler.get1d(); }
		}

		// write resulting images
		for (int iDim = 0; iDim < dimension; iDim++)
		{
			Fimage<double> f(imageSize);
			for (const Ivec2 & p : bound)
			{
				int index = Math::Index(bound, p);
				f.at(p) = currentMask[index * dimension + iDim];
			}
			FimageIo::Save(f, "white" + std::to_string(iDim) + ".pfm");
		}
		double currentEnergy = ComputeEnergy(currentMask, imageSize, dimension);

		int numMaxIterations = 100000;
		// simulate annealing
		for (int i = 0; i < numMaxIterations; i++)
		{
			if (i % 1000 == 0) std::cout << i << std::endl;

			// choose two pixels to randomly swap
			const Ivec2 p1 = Ivec2(sampler.get2d() * Vec2(imageSize));
			const Ivec2 p2 = Ivec2(sampler.get2d() * Vec2(imageSize));

			// randomly swap
			SwapPixel(&currentMask, imageSize, dimension, p1, p2);

			// compute the level of energy
			const double newEnergy = ComputeEnergy(currentMask, imageSize, dimension);

			// probablisiticly accept a new state
			if (newEnergy < currentEnergy)
			{
				currentEnergy = newEnergy;
			}
			else
			{
				// switch back to previous state
				SwapPixel(&currentMask, imageSize, dimension, p1, p2);
			}

			// write resulting images
			if (i % 1000 == 0 || i == numMaxIterations - 1)
				for (int iDim = 0; iDim < dimension; iDim++)
				{
					Fimage<double> f(imageSize);
					for (const Ivec2 & p : bound)
					{
						int index = Math::Index(bound, p);
						f.at(p) = currentMask[index * dimension + iDim];
					}
					FimageIo::Save(f, "blue" + std::to_string(i) + "_" + std::to_string(iDim) + ".pfm");
				}
		}
	}
};
