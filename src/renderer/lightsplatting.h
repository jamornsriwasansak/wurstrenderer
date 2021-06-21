#pragma once

#include "common/renderer.h"
#include "common/wurst.h"
#include "common/scene.h"
#include "common/sampler.h"
#include "common/camera.h"
#include "common/parallel.h"

#include "pathsampler/generalsubpathsampler.h"

struct LightSplattingRenderer : public SingleCameraRenderer
{
	static shared_ptr<LightSplattingRenderer> FromSimple2Json(shared_ptr<const Camera> & camera,
													   shared_ptr<const Scene> & scene,
													   shared_ptr<Sampler> & sampler,
													   const json & json,
													   const filesystem::path & rootJsonPath)
	{
		const int numVertices = Util::Json::GetValue<int>(json, "num_path_vertices", 10);
		const int numSpp = Util::Json::GetValue<int>(json, "num_spp", 16);
		return make_shared<LightSplattingRenderer>(camera, scene, sampler, numVertices, numSpp);
	}

	LightSplattingRenderer(shared_ptr<const Camera> & camera,
						   shared_ptr<const Scene> & scene,
						   shared_ptr<Sampler> & sampler,
						   const int numPathVertices,
						   const int numSamples) :
		SingleCameraRenderer(camera, scene, sampler),
		mNumPathVertices(numPathVertices),
		mNumSamples(numSamples),
		subpathSampler(mScene, false)
	{
	}

	void render() override
	{
		mCamera->mFilm->requestAtomicBuffer();

		// true number of samples = resolution[0] * resolution[1] * predefined numSamples
		double invTrueNumSamples = 1.0 / double(mNumSamples) / double(Math::Volume(mCamera->mFilm->mResolution));

		// create progress report and inform that PTNEE has progress report regarding this film
		shared_ptr<ProgressReport> progressReport = make_shared<ProgressReport>(Math::Volume(mCamera->mFilm->mResolution));
		Viewer::SetProgressReport(mCamera->mFilm, progressReport);

		Parallel::Split2d(mCamera->mFilm->mResolution, [&](const Ibound2 & bound)
		{
			int tileSeed = bound.pMin[0] * mCamera->mFilm->mResolution[1] + bound.pMin[1];
			shared_ptr<Sampler> tileSampler = mSampler->clone(tileSeed);

			for (const Ivec2 i : bound)
			{
				for (int i = 0; i < mNumSamples; i++)
				{
					if (mNumPathVertices < 2) continue;

					shared_ptr<CameraVertex> c0;

					// sample a camera vertex
					subpathSampler.createCameraSubpath(tileSampler.get(),
						nullptr,
						*mCamera,
						Ivec2(0, 0),
						1,
						SubpathSampler::Request::None,
						[&](shared_ptr<CameraVertex> c, const double) -> bool
						{
							c0 = c;
							return true;
						},
						nullptr);

					// trace light subpath
					subpathSampler.createLightSubpath(
						tileSampler.get(),
						mNumPathVertices - 1,
						SubpathSampler::Request::None,
						[&](shared_ptr<LightVertex> l, const double) -> bool
						{
							Vec2 raster;
							const double geometryTerm = Vertex::GeometryTerm(*l, *c0); if (geometryTerm == 0.0) return true;
							const Spectrum c0scatter = c0->scatter(&raster, CameraImaginaryVertex(), *l); if (Math::IsZero(c0scatter)) return true;
							const Spectrum lscatter = l->scatter(LightImaginaryVertex(), *c0); if (Math::IsZero(lscatter)) return true;
							const Spectrum transmittance = mScene->transmittance(*l, *c0); if (Math::IsZero(transmittance)) return true;

							const Spectrum r = c0->mPathContrib * c0scatter * geometryTerm * transmittance * lscatter * l->mPathContrib;
							mCamera->mFilm->atomicAddSample(raster, r * invTrueNumSamples);
							return true;
						},
						[&](shared_ptr<Vertex> l, shared_ptr<Vertex> lprev, const Vec3 &, const double) -> bool
						{
							if (l == nullptr) return false;
							if (l->hasType(VertexType::Surface | VertexType::Medium))
							{
								Vec2 raster;
								const double geometryTerm = Vertex::GeometryTerm(*l, *c0); if (geometryTerm == 0.0) return true;
								const Spectrum c0scatter = c0->scatter(&raster, CameraImaginaryVertex(), *l); if (Math::IsZero(c0scatter)) return true;
								const Spectrum lscatter = l->scatter(*lprev, *c0); if (Math::IsZero(lscatter)) return true;
								const Spectrum transmittance = mScene->transmittance(*l, *c0); if (Math::IsZero(transmittance)) return true;

								const Spectrum r = c0->mPathContrib * c0scatter * geometryTerm * transmittance * lscatter * l->mPathContrib;
								mCamera->mFilm->atomicAddSample(raster, r * invTrueNumSamples);
							}
							return true;
						});
				} // end num samples
				progressReport->increment(1);
			}
			Viewer::Redraw(mCamera->mFilm);
		}); // end parallel 2d
	}

	GeneralSubpathSampler subpathSampler;
	int mNumPathVertices;
	int mNumSamples;
};
