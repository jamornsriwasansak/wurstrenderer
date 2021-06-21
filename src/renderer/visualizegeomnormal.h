#pragma once

#include "common/renderer.h"
#include "common/wurst.h"
#include "common/scene.h"
#include "common/sampler.h"
#include "common/camera.h"
#include "common/parallel.h"

#include "pathsampler/generalsubpathsampler.h"

#include "geometry/plane.h"

struct VisualizeGeometryNormal : public SingleCameraRenderer
{
	static shared_ptr<VisualizeGeometryNormal> FromSimple2Json(shared_ptr<const Camera> & camera,
															   shared_ptr<const Scene> & scene,
															   shared_ptr<Sampler> & sampler,
															   const json & json,
															   const filesystem::path & rootJsonPath)
	{
		const int numSpp = Util::Json::GetValue<int>(json, "num_spp", 16);
		return make_shared<VisualizeGeometryNormal>(camera, scene, sampler, numSpp);
	}

	VisualizeGeometryNormal(shared_ptr<const Camera> & camera, shared_ptr<const Scene> & scene, shared_ptr<Sampler> & sampler, const int numSamples) :
		SingleCameraRenderer(camera, scene, sampler),
		mNumSamples(numSamples)
	{
	}

	void render() override
	{
		mCamera->mFilm->requestBaseBuffer();
		GeneralSubpathSampler pathSampler(mScene, true);
		double invNumSamples = 1.0 / static_cast<double>(mNumSamples);

		Parallel::Split2d(mCamera->mFilm->mResolution, [&](const Ibound2 & bound) {
			int tileSeed = bound.pMin[0] * mCamera->mFilm->mResolution[1] + bound.pMin[1];
			shared_ptr<Sampler> tileSampler = mSampler->clone(tileSeed);

			for (const Ivec2 pixel : bound)
			{
				Spectrum result(0.0);
				for (int iSample = 0; iSample < mNumSamples; iSample++)
				{
					pathSampler.createCameraSubpath(
						tileSampler.get(),
						nullptr,
						*mCamera,
						pixel,
						2,
						SubpathSampler::Request::None,
						nullptr,
						[&](shared_ptr<Vertex> & v, shared_ptr<Vertex> & prevVertex, const Vec3 &, const double) -> bool
						{
							if (v == nullptr) return false;
							if (v->hasType(VertexType::Light))
							{
								result += Spectrum(1.0);
							}
							else if (v->hasType(VertexType::Surface))
							{
								result += SpectrumConvertUtil::SpectrumFromRgb(v->mGeometryNormal);
							}
							return true;
						}
					);
				}
				mCamera->mFilm->addSample(pixel, result * invNumSamples);
			}
			Viewer::Redraw(mCamera->mFilm);
		});
	}

	Uint mNumSamples;
};
