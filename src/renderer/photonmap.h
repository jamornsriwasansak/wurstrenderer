#pragma once

#include "common/wurst.h"

#include "common/camera.h"
#include "common/hashgrid.h"
#include "common/parallel.h"
#include "common/renderer.h"
#include "common/scene.h"
#include "common/sampler.h"

#include "pathsampler/generalsubpathsampler.h"

struct PhotonMapRenderer : public SingleCameraRenderer
{
	static shared_ptr<PhotonMapRenderer> FromSimple2Json(shared_ptr<const Camera> & camera,
														 shared_ptr<const Scene> & scene,
														 shared_ptr<Sampler> & sampler,
														 const json & json,
														 const filesystem::path & rootJsonPath)
	{
		const int numPathVertices = Util::Json::GetValue<int>(json, "num_path_vertices", 10);
		const int numSpp = Util::Json::GetValue<int>(json, "num_spp", 16);
		const int numIterations = Util::Json::GetValue<int>(json, "num_iterations", 1);
		return make_shared<PhotonMapRenderer>(camera, scene, sampler, numPathVertices, numSpp);
	}

	PhotonMapRenderer(shared_ptr<const Camera> & camera,
					  shared_ptr<const Scene> & scene,
					  shared_ptr<Sampler> & sampler,
					  const int numPathVertices,
					  const int numSpp):
		SingleCameraRenderer(camera, scene, sampler),
		mNumPathVertices(numPathVertices),
		mNumLightSubpaths(numSpp * Math::Volume(camera->mFilm->mResolution)),
		mNumSpp(numSpp),
		mSubpathSampler(scene, true)
	{
	}

	struct Hitpoint
	{
		Hitpoint(const shared_ptr<Vertex> & vertex, const Ivec2 & pixel, const Vec3 & dirToCameraVertex): mVertex(vertex), mPixel(pixel), mDirToCameraVertex(dirToCameraVertex)
		{
		}

		shared_ptr<Vertex> mVertex;
		Ivec2 mPixel;
		Vec3 mDirToCameraVertex;
	};

	void render() override
	{
		double radius = 0.01;

		mCamera->mFilm->requestBaseBuffer();
		mCamera->mFilm->requestAtomicBuffer();

		int numPathsPerIteration = Math::Volume(mCamera->mFilm->mResolution);
		int numPhotonPathsPerIteration = 10000000;

		// create hitpoint grid
		Bound3 sceneBound = Bound3::Padding(mScene->bound(), radius * 2.0);
		Ivec3 resolution = Ivec3((sceneBound.pMax - sceneBound.pMin) / radius + Vec3(0.5));

		int photonRngOffset = Math::Volume(mCamera->mFilm->mResolution) * mIterations;

		for (int iTer = 0;iTer < mIterations;iTer++)
		{ 
			HashGrid<Hitpoint> hitpointGrid(sceneBound, resolution, numPathsPerIteration);

			// trace camera subpaths and store in hashgrid
			Parallel::Split2d(mCamera->mFilm->mResolution, [&](const Ibound2 & bound)
			{
				const int tileSeed = bound.pMin[1] * mCamera->mFilm->mResolution[0] + bound.pMin[0];
				const int offsetSeed = Math::Volume(mCamera->mFilm->mResolution) * iTer;
				shared_ptr<Sampler> tileSampler = mSampler->clone(tileSeed + offsetSeed);

				for (const Ivec2 pixel : bound)
				{
					// store camera hitpoint in grid
					Vec2 raster;
					mSubpathSampler.createCameraSubpath(tileSampler.get(), &raster, *mCamera, pixel,
						2, SubpathSampler::Request::None,
						nullptr,
						[&](shared_ptr<Vertex> c, shared_ptr<Vertex> cprev, const Vec3 &, const double) -> bool
						{
							if (c == nullptr) return false;
							if (c->hasType(VertexType::Surface))
							{
								hitpointGrid.insert(c->mPosition, radius, Hitpoint(c, pixel, Math::Normalize(*cprev - *c)));
								return true;
							}
							else if (c->hasType(VertexType::Light))
							{
								// add result
								const Spectrum le = c->mPathContrib * c->scatter(LightImaginaryVertex(), *cprev) * c->source();
								mCamera->mFilm->addSample(raster, le);
								return false;
							}
							return false;
						}
					);
				}
			});
			mCamera->mFilm->mBaseScalingFactor = 1.0 / double(iTer + 1);
			Viewer::Redraw(mCamera->mFilm);
			// end hitpoint creation

			// start neighbour searching
			Parallel::Split(0, numPhotonPathsPerIteration, [&](const int photonIndexStart, const int photonIndexEnd)
			{
				const int offsetSeed = numPhotonPathsPerIteration * iTer;
				shared_ptr<Sampler> photonSampler = mSampler->clone(photonRngOffset + photonIndexStart + offsetSeed);
				for (int iPhotonPath = photonIndexStart; iPhotonPath < photonIndexEnd; iPhotonPath++)
				{ 
					// trace photon
					mSubpathSampler.createLightSubpath(photonSampler.get(),
						mNumPathVertices - 1, SubpathSampler::Request::None,
						nullptr,
						[&](shared_ptr<Vertex> l, shared_ptr<Vertex> lprev, const Vec3 &, const double) -> bool
						{
							if (l == nullptr) return false;
							if (l->hasType(VertexType::Surface))
							{
								auto hashBucket = hitpointGrid.getSearchBucket(l->mPosition);
								for (const Hitpoint & hitpoint : hashBucket)
								{
									if (Math::Distance2(hitpoint.mVertex->mPosition, l->mPosition) <= radius * radius)
									{
										Spectrum scatter = hitpoint.mVertex->scatter(Math::Normalize(*lprev - *hitpoint.mVertex), hitpoint.mDirToCameraVertex);
										Spectrum contrib = hitpoint.mVertex->mPathContrib * scatter * l->mPathContrib;
										mCamera->mFilm->atomicAddSample(hitpoint.mPixel, contrib);
									}
								}
								// end neighbour search per photon
								return true;
							}
							return false;
						}
					);
				}
			});
			mCamera->mFilm->mAtomicScalingFactor = 1.0 / (Math::Pi * radius * radius * static_cast<double>(numPhotonPathsPerIteration * (iTer + 1)));
			Viewer::Redraw(mCamera->mFilm);
			// end neighbour searching
		}
	}

	GeneralSubpathSampler mSubpathSampler;
	int mIterations = 1;
	int mNumPathVertices;
	int mNumLightSubpaths;
	int mNumSpp;
};