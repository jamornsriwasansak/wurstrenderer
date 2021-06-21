#pragma once

#ifndef BDPT_MIS_WEIGHT_FAST
#define BDPT_MIS_WEIGHT_FAST
#endif

#include "common/renderer.h"
#include "common/wurst.h"
#include "common/scene.h"
#include "common/sampler.h"
#include "common/camera.h"
#include "common/parallel.h"
#include "common/viewer.h"
#include "common/scopedassign.h"

#include "pathsampler/generalsubpathsampler.h"

struct BdptSplat
{
	BdptSplat()
	{
	}

	BdptSplat(const Vec2 & raster, const Spectrum & contrib, const int s, const int t) : mRaster(raster), mContrib(contrib), mS(s), mT(t) {}

	Vec2 mRaster;
	Spectrum mContrib;
	int mS;
	int mT;
};

struct BdptVertex
{
	shared_ptr<Vertex> mVertex;
	double mPdfFromLight;
	double mPdfFromCamera;

	BdptVertex(const shared_ptr<Vertex> & vertex, const double pdfFromLight, const double pdfFromCamera):
		mVertex(vertex), mPdfFromLight(pdfFromLight), mPdfFromCamera(pdfFromCamera)
	{
	}
};

struct BidirPathRenderer : public SingleCameraRenderer
{
	static shared_ptr<BidirPathRenderer> FromSimple2Json(shared_ptr<const Camera> & camera,
														 shared_ptr<const Scene> & scene,
														 shared_ptr<Sampler> & sampler,
														 const json & json,
														 const filesystem::path & rootJsonPath)
	{
		const int numPathVertices = Util::Json::GetValue<int>(json, "num_path_vertices", 10);
		const int numSpp = Util::Json::GetValue<int>(json, "num_spp", 16);
		return make_shared<BidirPathRenderer>(camera, scene, sampler, numPathVertices, numSpp);
	}

	BidirPathRenderer(shared_ptr<const Camera> & camera,
					  shared_ptr<const Scene> & scene,
					  shared_ptr<Sampler> & sampler,
					  const int numVertices,
					  const int numSpp,
					  const bool doStratification = true) :
		SingleCameraRenderer(camera, scene, sampler),
		mNumPathVertices(numVertices),
		mNumSpp(numSpp),
		mSubpathSampler(scene, doStratification)
	{
	}

	int indexFromTechnique(const int s, const int t) const
	{
		assert(s >= 0);
		assert(t > 0);

		const int pathLength = s + t;
		const int offset = (pathLength * (pathLength - 1)) / 2 - 1;
		return offset + s;
	}

	Spectrum bdptConnect(Vec2 * raster,
						 const std::vector<BdptVertex> & lightInfos,
						 const std::vector<BdptVertex> & cameraInfos,
						 const int s,
						 const int t) const
	{
		assert(s >= 0);
		assert(t >= 1);
		assert(s + t >= 2);

		if (s == 0)
		{
			assert(t >= 2);
			Vertex & cp = *cameraInfos[t - 1].mVertex;
			Vertex & c = *cameraInfos[t].mVertex;

			if (!c.hasType(VertexType::Light)) return Spectrum(0.0);

			// compute splats
			return c.mPathContrib * c.scatter(LightImaginaryVertex(), cp) * c.source();
		}
		else
		{
			Vertex & lp = *lightInfos[s - 1].mVertex;
			Vertex & l = *lightInfos[s].mVertex;

			Vertex & cp = *cameraInfos[t - 1].mVertex;
			Vertex & c = *cameraInfos[t].mVertex;

			if (t == 1 && !c.hasType(VertexType::Camera)) return Spectrum(0.0);
			else if (t > 1 && !c.hasType(VertexType::Surface | VertexType::Medium)) return Spectrum(0.0);

			// compute splats
			const double geometryTerm = Vertex::GeometryTerm(l, c); if (geometryTerm == 0.0) return Spectrum(0.0);
			const Spectrum cscatter = c.scatter(raster, cp, l); if (Math::IsZero(cscatter)) return Spectrum(0.0);
			const Spectrum lscatter = l.scatter(lp, c); if (Math::IsZero(lscatter)) return Spectrum(0.0);
			const Spectrum transmittance = mScene->transmittance(l, c); if (Math::IsZero(transmittance)) return Spectrum(0.0);
			return l.mPathContrib * lscatter * transmittance * geometryTerm * cscatter * c.mPathContrib;
		}

		return Spectrum(0.0);
	}

	double computeMisWeight(std::vector<BdptVertex> & lightPoints,
							std::vector<BdptVertex> & cameraPoints,
							const int s,
							const int t,
							const int power) const
	{
		assert(s >= 0);
		assert(t >= 1);
		assert(s + t >= 2);

		const int numVertices = s + t;

		// suppose numVertices = 5, s = 3, t = 2, lightPoints indices: 0, 1, 2, 3, cameraPoints indices: 6, 5, 4
		auto getConnectedInfo = [&](const int iVertex) -> BdptVertex& { return (iVertex <= s) ? lightPoints[iVertex] : cameraPoints[numVertices - iVertex + 1]; };
		auto vertex = [&](const int iVertex) { return (iVertex <= s) ? lightPoints[iVertex].mVertex : cameraPoints[numVertices - iVertex + 1].mVertex; };
		auto lPdf = [&](const int iVertex) { return (iVertex <= s) ? lightPoints[iVertex].mPdfFromLight : cameraPoints[numVertices - iVertex + 1].mPdfFromLight; };
		auto cPdf = [&](const int iVertex) { return (iVertex <= s) ? lightPoints[iVertex].mPdfFromCamera : cameraPoints[numVertices - iVertex + 1].mPdfFromCamera; };

		// check whether the first vertex is on light source or not
		const shared_ptr<LightVertex> lv = Vertex::SafeCast<LightVertex>(getConnectedInfo(1).mVertex);
		if (!getConnectedInfo(1).mVertex->hasType(VertexType::Light)) return 0.0;

		ScopedAssignment<double> sa1, sa2, sa3, sa4;

		// compute pdf at the connection point
		BdptVertex 
			* c = &getConnectedInfo(s + 1), // last camera vertex
			* cp = &getConnectedInfo(s + 2), // penultimate camera vertex
			* l = (s >= 1) ? &getConnectedInfo(s) : nullptr, // last light vertex
			* lp = (s >= 1) ? &getConnectedInfo(s - 1) : nullptr; // penultimate light vertex

		// update prob
		if (s >= 1) { sa1 = { &l->mPdfFromCamera, c->mVertex->pdfA(*cp->mVertex, *l->mVertex) }; }
		if (s >= 2) { sa2 = { &lp->mPdfFromCamera, l->mVertex->pdfA(*c->mVertex, *lp->mVertex) }; }
		if (s >= 1 && t >= 2) { sa3 = { &c->mPdfFromLight, l->mVertex->pdfA(*lp->mVertex, *c->mVertex) }; }
		if (s >= 1 && t >= 2) { sa4 = { &cp->mPdfFromLight, c->mVertex->pdfA(*l->mVertex, *cp->mVertex) }; }

		auto computePathPdfA = [&](const int p, const int q) -> double
		{
			double lightPathPdfA = 1.0;
			double cameraPathPdfA = 1.0;
			for (int i = 1; i <= p; i++) lightPathPdfA *= lPdf(i);
			for (int i = p + 1; i < (p + q); i++) cameraPathPdfA *= cPdf(i);
			return lightPathPdfA * cameraPathPdfA;
		};

		double invMisWeight = 1.0;
		double lRatio = 1.0;

		for (int i = s; i > 0; i--)
		{
			#ifdef BDPT_MIS_WEIGHT_FAST
				lRatio *= cPdf(i) / lPdf(i);
			#else
				lRatio = computePathPdfA(i - 1, numVertices - (i - 1)) / computePathPdfA(s, t);
			#endif
			invMisWeight += std::pow(lRatio, power);
		}

		double cRatio = 1.0;
		for (int i = s + 1; i < numVertices; i++)
		{
			#ifdef BDPT_MIS_WEIGHT_FAST
				cRatio *= lPdf(i) / cPdf(i);
			#else
				cRatio = computePathPdfA(i, numVertices - i) / computePathPdfA(s, t);
			#endif
			invMisWeight += std::pow(cRatio, power);
		}

		double misWeight = 1.0 / invMisWeight;
		assert(misWeight <= 1.0);
		return misWeight;
	}

	void traceLightSubpath(std::vector<BdptVertex> * lVertices, Sampler * sampler) const
	{
		mSubpathSampler.createLightSubpath(
			sampler,
			mNumPathVertices,
			SubpathSampler::Request::PdfA,
			[&](shared_ptr<LightVertex> & vertex, const double pdfA) -> bool
			{ 
				lVertices->emplace_back(vertex, pdfA, 0.0);
				return true;
			},
			[&](shared_ptr<Vertex> & vertex, const shared_ptr<Vertex> & prevVertex, const Vec3 &, const double pdfA) -> bool
			{
				if (vertex == nullptr) return false;

				if (vertex->hasType(VertexType::Surface | VertexType::Medium))
				{
					if (lVertices->size() > 2)
					{
						// only pushback backward pdf if prev prev vertex is not lightimaginary vertex
						BdptVertex & prevPrevInfo = lVertices->at(lVertices->size() - 2);
						prevPrevInfo.mPdfFromCamera = prevVertex->pdfA(*vertex, *prevPrevInfo.mVertex);
					}
					lVertices->emplace_back(vertex, pdfA, 0.0);
				}
				return true;
			});
	}

	void traceCameraSubpath(std::vector<BdptVertex> * cameraVertices, Sampler * sampler, Vec2 * raster, const Ivec2 & pixel) const
	{
		// trace camera subpath and store in camera vertices vector
		mSubpathSampler.createCameraSubpath(
			sampler,
			raster,
			*mCamera,
			pixel,
			mNumPathVertices,
			SubpathSampler::Request::PdfA,
			[&](const shared_ptr<CameraVertex> & vertex, const double pdfA) -> bool
			{ 
				cameraVertices->emplace_back(vertex, pdfA, 0.0);
				return true;
			},
			[&](const shared_ptr<Vertex> & vertex, const shared_ptr<Vertex> & prevVertex, const Vec3 & dirToCurrentVertex, const double pdfA) -> bool
			{
				if (vertex == nullptr) return false;

				if (cameraVertices->size() > 2)
				{
					// only pushback backward pdf if prev prev vertex is not cameraimaginary vertex
					BdptVertex & prevPrevInfo = cameraVertices->at(cameraVertices->size() - 2);
					prevPrevInfo.mPdfFromLight = prevVertex->pdfA(*vertex, *prevPrevInfo.mVertex);
				}
				if (vertex->hasType(VertexType::Light))
				{
					// compute prevInfo bwdPdf
					BdptVertex & prevInfo = cameraVertices->at(cameraVertices->size() - 1);
					prevInfo.mPdfFromLight = vertex->pdfA(LightImaginaryVertex(), *prevInfo.mVertex);

					// compute bwdPdf of this vertex
					const double bwdPdfA = mSubpathSampler.evalLightPdfA(*Vertex::SafeCast<LightVertex>(vertex));
					cameraVertices->emplace_back(vertex, bwdPdfA, pdfA);
				}
				else if (vertex->hasType(VertexType::Surface | VertexType::Medium))
				{
					cameraVertices->emplace_back(vertex, 0.0, pdfA);
				}
				return true;
			});
	}

	void evalContribution(std::vector<BdptSplat> * splatContribs,
						  std::vector<BdptVertex> & lightPoints,
						  std::vector<BdptVertex> & cameraPoints,
						  const Vec2 & cameraRaster,
						  const int specifiedS = -1,
						  const int specifiedT = -1) const
	{
		int numS = static_cast<int>(lightPoints.size()) - 1;
		int numT = static_cast<int>(cameraPoints.size()) - 1;

		BdptSplat baseContrib(cameraRaster, Spectrum(0.0), 0, 0);

		for (int s = 0; s <= numS; s++)
			for (int t = 1; t <= numT; t++)
			{
				if (s + t <= 1 || s + t > mNumPathVertices) continue;
				if (specifiedS != -1 && specifiedS != s) continue;
				if (specifiedT != -1 && specifiedT != t) continue;

				// compute index for unweightedFilm
				int index = indexFromTechnique(s, t);

				// compute splats and splat on unweightedFilm
				Vec2 raster = cameraRaster;
				Spectrum result = bdptConnect(&raster, lightPoints, cameraPoints, s, t);
				if (Math::IsZero(result)) continue;

				double misWeight = computeMisWeight(lightPoints, cameraPoints, s, t, 1);

				if (t > 1)
				{
					if (mNeedSeperatedContrib || mDebugWeighted)
						splatContribs->emplace_back(raster, misWeight * result, s, t);
					else
						baseContrib.mContrib += misWeight * result;
				}
				else
				{
					// t == 1, special case. required for light splatting
					splatContribs->emplace_back(raster, misWeight * result, s, t);
				}
			}

		if (!(mNeedSeperatedContrib || mDebugWeighted)) splatContribs->emplace_back(baseContrib);
	}
	
	void render() override
	{
		if (mNumPathVertices <= 1) return;

		const double invNumSamples = 1.0 / double(mNumSpp);
		const double invNumLightSubpaths = 1.0 / double(Math::Volume(mCamera->mFilm->mResolution)) / double(mNumSpp);
		const int numTechniques = (mNumPathVertices * (mNumPathVertices + 1)) / 2 - 1;
		mCamera->mFilm->mBaseScalingFactor = invNumSamples;
		mCamera->mFilm->mAtomicScalingFactor = invNumLightSubpaths;
		mCamera->mFilm->requestAtomicBuffer();
		mCamera->mFilm->requestBaseBuffer();

		// create progress report
		shared_ptr<ProgressReport> progressReport = make_shared<ProgressReport>(mNumSpp, Math::Volume(mCamera->mFilm->mResolution));
		Viewer::SetProgressReport(mCamera->mFilm, progressReport);

		// debugging films for different techniques
		if (mDebugWeighted)
			for (int pathLength = 2; pathLength <= mNumPathVertices; pathLength++)
				for (int s = 0; s < pathLength; s++)
				{
					shared_ptr<Film> weightedFilm = make_shared<Film>(mCamera->mFilm->mResolution);
					weightedFilm->mBaseScalingFactor = invNumSamples;
					weightedFilm->mAtomicScalingFactor = invNumLightSubpaths;
					if (pathLength - s == 1)
						weightedFilm->requestAtomicBuffer();
					else
						weightedFilm->requestBaseBuffer();
					mWeightedFilms.push_back(weightedFilm);
				}

		std::atomic<int> numDrawnTiles = 0;

		Parallel::Split2d(mCamera->mFilm->mResolution, [&](const Ibound2 & bound)
		{
			const int tileSeed = bound.pMin[0] * mCamera->mFilm->mResolution[1] + bound.pMin[1];
			shared_ptr<Sampler> tileSampler = mSampler->clone(tileSeed);

			// create vector for light subpath vertices
			std::vector<BdptVertex> lVertices;
			lVertices.reserve(mNumPathVertices + 1);

			// create vector for camera subpath vertices
			std::vector<BdptVertex> cVertices;
			cVertices.reserve(mNumPathVertices + 1);

			// contributions of all techniques
			std::vector<BdptSplat> splats(numTechniques, BdptSplat());

			for (const Ivec2 & pixel : bound)
			{
				for (int iSample = 0; iSample < mNumSpp; iSample++)
				{
					// trace light subpath
					lVertices.clear();
					lVertices.emplace_back(make_shared<LightImaginaryVertex>(), 0.0, 0.0);
					traceLightSubpath(&lVertices, tileSampler.get());

					// start preparing vector for tracing from camera
					Vec2 cameraRaster;
					cVertices.clear();
					cVertices.emplace_back(make_shared<CameraImaginaryVertex>(), 0.0, 0.0);
					traceCameraSubpath(&cVertices, tileSampler.get(), &cameraRaster, pixel);

					// eval and add to film
					splats.clear();
					evalContribution(&splats, lVertices, cVertices, cameraRaster);
					for (BdptSplat & c : splats)
					{
						if (c.mT == 1)
							mCamera->mFilm->atomicAddSample(c.mRaster, c.mContrib);
						else
							mCamera->mFilm->addSample(c.mRaster, c.mContrib);

						if (mDebugWeighted)
						{
							int index = indexFromTechnique(c.mS, c.mT);
							if (c.mT == 1)
								mWeightedFilms[index]->atomicAddSample(c.mRaster, c.mContrib);
							else
								mWeightedFilms[index]->addSample(c.mRaster, c.mContrib);
						}
					}
				}

				// report progress to viewer
				progressReport->increment(mNumSpp);
			};

			if (numDrawnTiles++ % 128 == 127)
				Viewer::Redraw(mCamera->mFilm);
			else
				Viewer::Redraw(mCamera->mFilm, bound);
		});

		// write debug film
		if (mDebugWeighted)
			for (int pathLength = 2; pathLength <= mNumPathVertices; pathLength++)
				for (int s = 0; s < pathLength; s++)
				{
					int t = pathLength - s;
					int index = indexFromTechnique(s, t);
					FimageIo::Save(mWeightedFilms[index]->getImage(), "debug_weighted_" + std::to_string(s + t) + "_" + std::to_string(s) + "_" + std::to_string(t) + ".pfm");
					FimageIo::Save(mWeightedFilms[index]->getImage(), "debug_weighted_" + std::to_string(s + t) + "_" + std::to_string(s) + "_" + std::to_string(t) + ".png");
				}
	}

	mutable std::vector<shared_ptr<Film>> mWeightedFilms;

	GeneralSubpathSampler mSubpathSampler;
	int mNumPathVertices;
	int mNumSpp;
	bool mDebugWeighted = false;
	bool mNeedSeperatedContrib = false;
};