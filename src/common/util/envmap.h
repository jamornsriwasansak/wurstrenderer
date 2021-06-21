#pragma once

#include "common/wurst.h"
#include "common/sampler.h"
#include "common/sphericalrepresentation.h"
#include "sphericalrepresentation/equirectangular.h"

namespace Util
{
	struct Envmap
	{
		static EquiRectRep<double> ComputeDiffuse(const SphRep<double> & sphrep, const Ivec2 & size, const int numSamples, Sampler * sampler)
		{
			assert(numSamples >= 0);
			Fimage<double> result(size);
			result.forEachPixel([&](const Vec2 & uv, double * color)
								{
									const CoordFrame3 basis(Mapping::WorldFromPanorama(uv));

									// cosine sampling
									//for (int iSample = 0; iSample < numSamples; iSample++)
									std::atomic<double> sum = 0.0;
									Parallel::Split(0, numSamples, [&](int start, int end)
													{
														for (int i = start; i < end; i++)
														{
															const Vec3 sample = Mapping::CosineWeightedHemisphereFromSquare(sampler->get2d());
															const Vec3 rotated = basis.toWorld(sample);
															double s = 0;
															double v = sphrep.eval(rotated);
															while (!sum.compare_exchange_weak(s, s + v));
														}
													});

									*color = sum / static_cast<double>(numSamples);
								});

			return result;
		}
	};
}
