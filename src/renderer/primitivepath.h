#pragma once

#include "common/wurst.h"

#include "common/camera.h"
#include "common/parallel.h"
#include "common/renderer.h"
#include "common/scene.h"
#include "common/sampler.h"
#include "common/util/json.h"

#include "pathsampler/generalsubpathsampler.h"

// a path tracing at its simplest form. only bounce materials around.
struct PrimitivePathRenderer : public SingleCameraRenderer
{
	static shared_ptr<PrimitivePathRenderer> FromSimple2Json(shared_ptr<const Camera> & camera,
															 shared_ptr<const Scene> & scene,
															 shared_ptr<Sampler> & sampler,
															 const json & json,
															 const filesystem::path & rootJsonPath)
	{
		const int numPathVertices = Util::Json::GetValue<int>(json, "num_path_vertices", 10);
		const int numSpp = Util::Json::GetValue<int>(json, "num_spp", 16);
		return make_shared<PrimitivePathRenderer>(camera, scene, sampler, numPathVertices, numSpp);
	}

	PrimitivePathRenderer(shared_ptr<const Camera> & camera,
						  shared_ptr<const Scene> & scene,
						  shared_ptr<Sampler> & sampler,
						  const int numPathVertices,
						  const int numSpp,
						  const bool doStratification = true):
		SingleCameraRenderer(camera, scene, sampler),
		mDoStratification(doStratification),
		mNumPathVertices(numPathVertices),
		mNumSpp(numSpp),
		mSubpathSampler(scene, doStratification, true)
	{
	}

	void evalContribution(Splat * splat, Sampler * sampler, const Ivec2 & pixel = Ivec2(0, 0)) const
	{
		bool passed = false;
		// trace camera subpath
		mSubpathSampler.createCameraSubpath(sampler, &splat->mRaster, *mCamera, pixel, mNumPathVertices, SubpathSampler::Request::None,
			nullptr,
			[&](shared_ptr<Vertex> & vertex, shared_ptr<Vertex> & prevVertex, const Vec3 & dirToPrevVertex, const double) -> bool
			{
				if (vertex == nullptr) return false;
				if (vertex->hasType(VertexType::Light))
				{
					splat->mContrib = vertex->mPathContrib * vertex->scatter(LightImaginaryVertex(), *prevVertex) * vertex->source();
					return false;
				}
				return true;
			});
	}

	void render() override
	{
		if (mDoStratification)
			mCamera->mFilm->requestBaseBuffer();
		else
			mCamera->mFilm->requestAtomicBuffer();

		shared_ptr<ProgressReport> progressReport = make_shared<ProgressReport>(mNumSpp, Math::Volume(mCamera->mFilm->mResolution));
		Viewer::SetProgressReport(mCamera->mFilm, progressReport);

		mCamera->mFilm->mBaseScalingFactor = 1.0 / double(mNumSpp);
		mCamera->mFilm->mAtomicScalingFactor = 1.0 / double(mNumSpp) / double(Math::Volume(mCamera->mFilm->mResolution));

		Parallel::Split2d(mCamera->mFilm->mResolution, [&](const Ibound2 & bound)
		{
			// create sampler
			int tileSeed = bound.pMin[0] * mCamera->mFilm->mResolution[1] + bound.pMin[1];
			shared_ptr<Sampler> tileSampler = mSampler->clone(tileSeed);

			// create contribution
			Spectrum result(0.0);
			for (const Ivec2 pixel : bound)
			{
				for (int iSample = 0; iSample < mNumSpp; iSample++)
				{
					Splat contrib;
					evalContribution(&contrib, tileSampler.get(), pixel);
					if (mDoStratification)
						mCamera->mFilm->addSample(contrib.mRaster, contrib.mContrib);
					else
						mCamera->mFilm->atomicAddSample(contrib.mRaster, contrib.mContrib);
				}
				progressReport->increment(mNumSpp);
			}
			Viewer::Redraw(mCamera->mFilm, bound);
		});
	}

	GeneralSubpathSampler mSubpathSampler;
	bool mDoStratification;
	int mNumPathVertices;
	int mNumSpp;
};
