#pragma once

#include "common/wurst.h"

#include "common/renderer.h"
#include "common/scene.h"
#include "common/sampler.h"
#include "common/camera.h"
#include "common/parallel.h"
#include "common/viewer.h"

#include "renderer/bidirpath.h"

struct MbdptLSubpath
{
	MbdptLSubpath()
	{
	}

	std::vector<BdptVertex> mVertices;
	double mSortWeight;
};

struct MbdptCSubpath
{
	MbdptCSubpath()
	{
	}

	std::vector<BdptVertex> mVertices;
	Vec2 mRaster;
	double mSortWeight;
};

struct MatrixBidirPathRenderer : public SingleCameraRenderer
{
	MatrixBidirPathRenderer(shared_ptr<const Camera> & camera,
							shared_ptr<const Scene> & scene,
							shared_ptr<Sampler> & sampler,
							const int numVertices,
							const int numSpp,
							const bool doStratification = true):
		SingleCameraRenderer(camera, scene, sampler),
		mNumPathVertices(numVertices),
		mNumSpp(numSpp),
		mBdpt(camera, scene, sampler, numVertices, numSpp, doStratification)
	{
		mBdpt.mDebugWeighted = mDebugWeighted;
	}

	static double ComputeSubpathLength(const std::vector<BdptVertex> & vertices, const int numST)
	{
		// return infinite if path is terminate before specified numST
		if (numST > static_cast<int>(vertices.size())) { return std::numeric_limits<double>::infinity(); }

		// else compute path length
		double length = 0.0;
		for (int i = 2; i < vertices.size(); i++)
		{
			length += Math::Distance(vertices[i - 1].mVertex->mPosition, vertices[i].mVertex->mPosition);
		}
		return length;
	}

	void render() override
	{
		if (mNumPathVertices <= 1) return;

		const double invNumSamples = 1.0 / static_cast<double>(mNumSpp);
		const int numTechniques = (mNumPathVertices * (mNumPathVertices + 1)) / 2 - 1;
		const int numPixels = Math::Volume(mCamera->mFilm->mResolution);
		mCamera->mFilm->mAtomicScalingFactor = invNumSamples;
		mCamera->mFilm->requestAtomicBuffer();

		// create progress report
		shared_ptr<ProgressReport> progressReport = make_shared<ProgressReport>(mNumSpp, Math::Volume(mCamera->mFilm->mResolution));
		Viewer::SetProgressReport(mCamera->mFilm, progressReport);

		// debugging films for different techniques
		if (mDebugWeighted)
			for (int pathLength = 2; pathLength <= mNumPathVertices; pathLength++)
				for (int s = 0; s < pathLength; s++)
				{
					shared_ptr<Film> weightedFilm = make_shared<Film>(mCamera->mFilm->mResolution);
					weightedFilm->mAtomicScalingFactor = invNumSamples;
					if (pathLength - s == 1)
					{
						weightedFilm->requestAtomicBuffer();
					}
					else
					{
						weightedFilm->requestBaseBuffer();
					}
					mWeightedFilms.push_back(weightedFilm);
				}

		std::atomic<int> numDrawnTiles = 0;

		for (int iTer = 0; iTer < mNumSpp; iTer++)
		{
			Parallel::Split2d(mCamera->mFilm->mResolution, [&](const Ibound2 & bound)
			{
				const int tileSeed = (bound.pMin[0] * mCamera->mFilm->mResolution[1] + bound.pMin[1]) + (iTer * numPixels);
				shared_ptr<Sampler> tileSampler = mSampler->clone(tileSeed);
				const int numPixelsPerTile = Ibound2::Volume(bound);

				// generate eye and light subpaths
				std::vector<MbdptLSubpath> lSubpaths(numPixelsPerTile);
				std::vector<MbdptCSubpath> cSubpaths(numPixelsPerTile);
				for (const Ivec2 & pixel : bound)
				{
					// trace light subpath
					MbdptLSubpath & lSubpath = lSubpaths[Math::Index(bound, pixel)];
					lSubpath.mVertices.clear();
					lSubpath.mVertices.emplace_back(make_shared<LightImaginaryVertex>(), 0.0, 0.0);
					mBdpt.traceLightSubpath(&lSubpath.mVertices, tileSampler.get());

					// trace camera subpath
					MbdptCSubpath & cSubpath = cSubpaths[Math::Index(bound, pixel)];
					cSubpath.mVertices.clear();
					cSubpath.mVertices.emplace_back(make_shared<CameraImaginaryVertex>(), 0.0, 0.0);
					mBdpt.traceCameraSubpath(&cSubpath.mVertices, tileSampler.get(), &cSubpath.mRaster, pixel);
				};

				// create connections
				// TODO:: replace this random sampler with low discreprancy sampler
				RandomSampler connectionSampler(tileSeed + iTer * 132851);
				std::vector<Ivec2> connections(numPixelsPerTile);
				for (int i = 0; i < numPixelsPerTile; i++) { connections[i] = Ivec2(connectionSampler.get2d() * static_cast<double>(numPixelsPerTile)); }

				std::vector<BdptSplat> bdptSplatTmp(1);
				for (int s = 0; s <= mNumPathVertices; s++)
					for (int t = 1; t <= mNumPathVertices; t++)
					{
						if (s + t <= 1 || s + t > mNumPathVertices) continue;

						// assign sort weights and sort eye and light vertices
						for (int i = 0; i < numPixelsPerTile; i++)
						{
							lSubpaths[i].mSortWeight = ComputeSubpathLength(lSubpaths[i].mVertices, s);
							cSubpaths[i].mSortWeight = ComputeSubpathLength(cSubpaths[i].mVertices, t);
						}
						std::sort(lSubpaths.begin(), lSubpaths.end(), [&](const MbdptLSubpath & a, const MbdptLSubpath & b) { return a.mSortWeight < b.mSortWeight; });
						std::sort(cSubpaths.begin(), cSubpaths.end(), [&](const MbdptCSubpath & a, const MbdptCSubpath & b) { return a.mSortWeight < b.mSortWeight; });

						double scale = 1.0;
						if (t == 1) scale = 1.0 / static_cast<double>(Math::Volume(mCamera->mFilm->mResolution));

						for (int i = 0; i < numPixelsPerTile; i++)
						{
							bdptSplatTmp.clear();
							MbdptLSubpath & lSubpath = lSubpaths[connections[i][0]];
							MbdptCSubpath & cSubpath = cSubpaths[connections[i][1]];
							mBdpt.evalContribution(&bdptSplatTmp, lSubpath.mVertices, cSubpath.mVertices, cSubpath.mRaster, s, t);
							assert(bdptSplatTmp.size() <= 1);
							if (bdptSplatTmp.size() == 1)
								mCamera->mFilm->atomicAddSample(bdptSplatTmp[0].mRaster, bdptSplatTmp[0].mContrib * scale);
						};
					}
				Viewer::Redraw(mCamera->mFilm, bound);
				progressReport->increment(numPixelsPerTile);
			});
		}

		// write debug film
		if (mDebugWeighted)
			for (int pathLength = 2; pathLength <= mNumPathVertices; pathLength++)
				for (int s = 0; s < pathLength; s++)
				{
					int t = pathLength - s;
					int index = mBdpt.indexFromTechnique(s, t);
					FimageIo::Save(mWeightedFilms[index]->getImage(), "debug_weighted_" + std::to_string(s + t) + "_" + std::to_string(s) + "_" + std::to_string(t) + ".pfm");
					FimageIo::Save(mWeightedFilms[index]->getImage(), "debug_weighted_" + std::to_string(s + t) + "_" + std::to_string(s) + "_" + std::to_string(t) + ".png");
				}
	}

	mutable std::vector<shared_ptr<Film>> mWeightedFilms;

	BidirPathRenderer mBdpt;
	bool mDebugWeighted = true;
	int mNumPathVertices;
	int mNumSpp;
};
