#pragma once

#include "common/renderer.h"
#include "common/wurst.h"
#include "common/scene.h"
#include "common/sampler.h"
#include "common/camera.h"
#include "common/parallel.h"
#include "common/viewer.h"

#include "pathsampler/generalsubpathsampler.h"

#define NEE_USE_MIS

struct PathNeeRenderer : public SingleCameraRenderer
{
	static shared_ptr<PathNeeRenderer> FromSimple2Json(shared_ptr<const Camera> & camera,
													   shared_ptr<const Scene> & scene,
													   shared_ptr<Sampler> & sampler,
													   const json & json,
													   const filesystem::path & rootJsonPath)
	{
		const int numVertices = Util::Json::GetValue<int>(json, "num_path_vertices", 10);
		const int numSpp = Util::Json::GetValue<int>(json, "num_spp", 16);
		return make_shared<PathNeeRenderer>(camera, scene, sampler, numVertices, numSpp);
	}

	PathNeeRenderer(shared_ptr<const Camera> & camera,
					shared_ptr<const Scene> & scene,
					shared_ptr<Sampler> & sampler,
					const int numPathVertices,
					const int numSpp,
					const bool doStratification = true) :
		SingleCameraRenderer(camera, scene, sampler),
		mNumPathVertices(numPathVertices),
		mNumSpp(numSpp),
		mSubpathSampler(scene, doStratification, true)
	{
	}

	void evalContribution(Splat * splats, Sampler * sampler, const Ivec2 & pixel = Ivec2(0, 0)) const
	{
		// trace light subpath
		shared_ptr<Vertex> l;
		double lPdfA = 1.0;

        #ifdef NEE_USE_MIS
            const SubpathSampler::Request request = SubpathSampler::Request::PdfA;
        #else
            const SubpathSampler::Request request = SubpathSampler::Request::None;
        #endif

		mSubpathSampler.createLightSubpath(
			sampler,
			1,
			request,
			[&](shared_ptr<LightVertex> & vertex, const double pdfA) -> bool
			{ 
				l = vertex;
				lPdfA = pdfA;
				return true;
			},
			nullptr);

		// trace camera subpath
		Uint cNumVertices = 1;

		mSubpathSampler.createCameraSubpath(
			sampler,
			&splats[0].mRaster,
			*mCamera,
			pixel,
			mNumPathVertices,
			request,
			[&](shared_ptr<CameraVertex> & c0, const double) -> bool
			{
				// compute splats
				Vec2 raster;
				const double geometryTerm    = Vertex::GeometryTerm(*l, *c0); if (geometryTerm == 0.0) return true;
				const Spectrum c0scatter     = c0->scatter(&raster, CameraImaginaryVertex(), *l); if (Math::IsZero(c0scatter)) return true;
				const Spectrum lscatter      = l->scatter(LightImaginaryVertex(), *c0); if (Math::IsZero(lscatter)) return true;
				const Spectrum transmittance = mScene->transmittance(*l, *c0); if (Math::IsZero(transmittance)) return true;

				// compute mis weight
				const double lightPdfA = lPdfA;
				const double camPdfA = c0->pdfA(CameraImaginaryVertex(), *l);
				#ifdef NEE_USE_MIS
					const double misWeight = lightPdfA / (lightPdfA + camPdfA);
				#else
					const double misWeight = 1.0;
				#endif
				assert(misWeight <= 1.0);
				const Spectrum r = misWeight * c0->mPathContrib * c0scatter * geometryTerm * transmittance * lscatter * l->mPathContrib;

				// light splatting splats
				splats[1].mRaster = raster;
				splats[1].mContrib = r;
				return true;
			},
			[&](shared_ptr<Vertex> & vertex, shared_ptr<Vertex> & prevVertex, const Vec3 &, const double cPdfA) -> bool
			{
				if (vertex == nullptr) return false;

				cNumVertices += 1;

				if (vertex->hasType(VertexType::Surface | VertexType::Medium) && (cNumVertices != mNumPathVertices))
				{
					// compute splats
					const double geometryTerm = Vertex::GeometryTerm(*vertex, *l); if (geometryTerm == 0.0) return true;
					const Spectrum transmittance = mScene->transmittance(*vertex, *l); if (Math::IsZero(transmittance)) return true;
					const Spectrum scatter1 = vertex->scatter(*prevVertex, *l); if (Math::IsZero(scatter1)) return true;
					const Spectrum scatter2 = l->scatter(LightImaginaryVertex(), *vertex); if (Math::IsZero(scatter2)) return true;

					// compute mis weight
					const double lightPdfA = lPdfA;
					const double bsdfPdfA = vertex->pdfA(*prevVertex, *l);
					#ifdef NEE_USE_MIS
						const double misWeight = lightPdfA / (lightPdfA + bsdfPdfA);
					#else
						const double misWeight = 1.0;
					#endif
					assert(misWeight <= 1.0);

					// add result
					splats[0].mContrib += misWeight * vertex->mPathContrib * scatter1 * scatter2 * transmittance * geometryTerm * l->mPathContrib;
					return true;
				}
				else if (vertex->hasType(VertexType::Light))
				{
					// compute mis weight
					const double bsdfPdfA = cPdfA; // for the second vertex, this is actually camera pdf
					const double lightPdfA = mSubpathSampler.evalLightPdfA(*Vertex::SafeCast<LightVertex>(vertex));
					#ifdef NEE_USE_MIS
						const double misWeight = bsdfPdfA / (lightPdfA + bsdfPdfA);
					#else
						const double misWeight = 0.0;
					#endif
					assert(misWeight <= 1.0);

					// add result
					splats[0].mContrib += misWeight
						* vertex->mPathContrib
						* vertex->scatter(LightImaginaryVertex(), *prevVertex)
						* vertex->source();
					return true;
				}

				return false;
			});
	}

	void render() override
	{
		mCamera->mFilm->requestAtomicBuffer();
		mCamera->mFilm->requestBaseBuffer();
		mCamera->mFilm->mBaseScalingFactor = 1.0 / double(mNumSpp);
		mCamera->mFilm->mAtomicScalingFactor = 1.0 / double(mNumSpp) / double(Math::Volume(mCamera->mFilm->mResolution));

		// create progress report and inform that PTNEE has progress report regarding this film
		shared_ptr<ProgressReport> progressReport = make_shared<ProgressReport>(mNumSpp, Math::Volume(mCamera->mFilm->mResolution));
		Viewer::SetProgressReport(mCamera->mFilm, progressReport);

		if (mNumPathVertices <= 1) return;

		std::atomic<int> numDrawnTiles = 0;

		Parallel::Split2d(mCamera->mFilm->mResolution, [&](const Ibound2 & bound)
		{
			const int tileSeed = bound.pMin[0] * mCamera->mFilm->mResolution[1] + bound.pMin[1];
			shared_ptr<Sampler> tileSampler = mSampler->clone(tileSeed);

			for (const Ivec2 pixel : bound)
			{
				tileSampler->setPixel(pixel);
				for (int iSample = 0; iSample < mNumSpp; iSample++)
				{
					Splat splats[2];
					evalContribution(splats, tileSampler.get(), pixel);
					mCamera->mFilm->addSample(splats[0].mRaster, splats[0].mContrib);
					mCamera->mFilm->atomicAddSample(splats[1].mRaster, splats[1].mContrib);
					tileSampler->nextSample();
				}
				progressReport->increment(mNumSpp);
			}

			if (numDrawnTiles++ % 128 == 127)
				Viewer::Redraw(mCamera->mFilm);
			else
				Viewer::Redraw(mCamera->mFilm, bound);
		});
	}

	GeneralSubpathSampler mSubpathSampler;
	int mNumPathVertices;
	int mNumSpp;
};
