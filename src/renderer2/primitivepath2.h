#pragma once

#include "common/wurst.h"

#include "common/camera2.h"
#include "common/parallel.h"
#include "common/renderer.h"
#include "common/scene.h"
#include "common/sampler.h"
#include "common/viewer.h"
#include "common/visualizer2.h"

#include "pathsampler2/simplesubpathsampler2.h"

// a path tracing at its simplest form. only bounce materials around.
struct PrimitivePathRenderer2 : public SingleCameraRenderer2
{
	static shared_ptr<PrimitivePathRenderer2> FromSimple2Json(shared_ptr<const Camera2> & camera,
															  shared_ptr<const Scene> & scene,
															  shared_ptr<Sampler> & sampler,
															  const json & json,
															  const filesystem::path & rootJsonPath)
	{
		const int numPathVertices = Util::Json::GetValue<int>(json, "num_path_vertices", 10);
		const int numSpp = Util::Json::GetValue<int>(json, "num_spp", 16);
		return make_shared<PrimitivePathRenderer2>(camera, scene, sampler, numPathVertices, numSpp, false);
	}

	PrimitivePathRenderer2(shared_ptr<const Camera2> & camera,
						   shared_ptr<const Scene> & scene,
						   shared_ptr<Sampler> & sampler,
						   const int numVertices,
						   const int numSpp,
						   const bool doVisualize) :
		SingleCameraRenderer2(camera, scene, sampler),
		mNumVertices(numVertices),
		mNumSpp(numSpp),
		mDoVisualize(doVisualize)
	{
		mCamera2->mFilm->request(FilmRequest::BaseBuffer);
	}

	void render() override
	{
		SimpleSubpathSampler2 subpathSampler(mScene, true);
		Visualizer2 pathVisualizer(mCamera2, mScene, Ivec2(512, 512));
		if (mDoVisualize) Viewer::Begin(pathVisualizer.mFilm);
		Parallel::Split2d(mCamera2->mFilm->mResolution, [&](const Ibound2 & bound)
		{
			int tileSeed = bound.pMin[0] * mCamera2->mFilm->mResolution[1] + bound.pMin[1];
			shared_ptr<Sampler> tileSampler = mSampler->clone(tileSeed);

			for (const Ivec2 pixel : bound)
			{
				Spectrum sum(0.0);
				for (Uint iSpp = 0; iSpp < mNumSpp; iSpp++)
				{
					Spectrum cContrib(1.0);
					Spectrum v(0.0);
					std::vector<shared_ptr<Vertex>> path;
					subpathSampler.createCameraSubpath(
						mSampler.get(),
						nullptr,
						*mCamera2,
						pixel,
						mNumVertices,
						[&](shared_ptr<CameraVertex> vertex, const Spectrum & contrib, const double)
						{
							cContrib *= contrib;
							if (mDoVisualize) path.push_back(vertex);
						},
						[&](shared_ptr<Vertex> vertex, shared_ptr<Vertex> prevVertex, const Spectrum & contrib, const double)
						{
							cContrib *= contrib;
							if (vertex->hasType(VertexType::Light))
							{
								//v = cContrib * vertex->scatter2(*prevVertex, LightImaginaryVertex()) * vertex->source2();
								v = cContrib;
							}
							if (mDoVisualize) path.push_back(vertex);
						}
					);
					if (mDoVisualize) pathVisualizer.drawPath(path);
					sum += v;
				}
				mCamera2->mFilm->addSample(pixel, sum / static_cast<double>(mNumSpp));
			}
			Viewer::Redraw(mCamera2->mFilm);
			if (mDoVisualize) Viewer::Redraw(pathVisualizer.mFilm);
		});
		if (mDoVisualize)
		{
			Viewer::End(pathVisualizer.mFilm);
			FimageIo::Save(pathVisualizer.mFilm->getImage(), "visualized.pfm");
		}
	}

	int mNumVertices;
	int mNumSpp;
	bool mDoVisualize;
};
