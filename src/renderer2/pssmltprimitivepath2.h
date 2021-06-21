#pragma once

#include "common/wurst.h"

#include "common/camera2.h"
#include "common/mlt.h"
#include "common/parallel.h"
#include "common/renderer.h"
#include "common/scene.h"
#include "common/sampler.h"
#include "common/visualizer2.h"
#include "common/viewer.h"

#include "pathsampler2/simplesubpathsampler2.h"

// a pssmlt over a primitive path tracer.

struct PssmltPrimitivePathRenderer2 : public SingleCameraRenderer2
{
	PssmltPrimitivePathRenderer2(shared_ptr<const Camera2> & camera,
							 shared_ptr<const Scene> & scene,
							 shared_ptr<Sampler> & sampler,
							 const Uint numVertices,
							 const Uint numSpp):
		SingleCameraRenderer2(camera, scene, sampler),
		mNumVertices(numVertices),
		mNumMutations(numSpp),
		mSubpathSampler2(scene, false)
	{
	}

	// basically doing primitive path tracing
	void evalContribution(PssmltState * state, std::vector<shared_ptr<Vertex>> * path) const
	{
		PssmltSampler replayableSampler(&state->mPrimarySample);

		state->mContrib = Spectrum(0.0);

		Spectrum cContrib(1.0);
		Vec2 pixel;

		mSubpathSampler2.createCameraSubpath(
			&replayableSampler,
			&state->mRaster,
			*mCamera2,
			Ivec2(0, 0),
			mNumVertices,
			[&](shared_ptr<CameraVertex> vertex, const Spectrum & c, const double)
			{
				cContrib *= c;
				if (path) path->push_back(vertex);
			},
			[&](shared_ptr<Vertex> vertex, shared_ptr<Vertex> prevVertex, const Spectrum & c, const double)
			{
				cContrib *= c;
				if (vertex->hasType(VertexType::Light))
				{
					state->mContrib = cContrib * vertex->scatter2(*prevVertex, LightImaginaryVertex()) * vertex->source2();
				}
				if (path) path->push_back(vertex);
			}
		);
		state->mScalarContrib = Math::Length(state->mContrib);
	}

	PssmltState perturbLarge(Rng * rng, std::vector<shared_ptr<Vertex>> * path) const
	{
		PssmltState result(mNumVertices);
		for (Uint i = 0; i < mNumVertices; i++)
		{
			result.mPrimarySample[i] = rng->nextFloat();
		}
		evalContribution(&result, path);
		return result;
	}

	PssmltState perturbSmall(const PssmltState & previous, Rng * rng, std::vector<shared_ptr<Vertex>> * path) const
	{
		PssmltState result(mNumVertices);
		for (Uint i = 0; i < mNumVertices; i++)
		{
			result.mPrimarySample[i] = MltUtil::PerturbExpLog(previous.mPrimarySample[i], 1.0 / 1024.0, 1.0 / 64.0, rng->nextFloat());
		}
		evalContribution(&result, path);
		return result;
	}

	void render() override
	{
		Visualizer2 visualizer(mCamera2, mScene, Ivec2(512, 512));
		Rng rng(0);

		std::vector<shared_ptr<Vertex>> path;
		std::vector<shared_ptr<Vertex>> * pathPtr = nullptr;
		if (mDoVisualize)
		{
			Viewer::Begin(visualizer.mFilm);
			pathPtr = &path;
		}

		// estimate normalization factor
		double b = 0.0;
		for (Uint iB = 0; iB < 100000; iB++)
		{
			b += perturbLarge(&rng, nullptr).mScalarContrib;
		}
		b /= 100000.0;

		PssmltState current = perturbLarge(&rng, pathPtr);
		for (Uint iMutation = 0; iMutation < mNumMutations; iMutation++)
		{
			PssmltState proposed = (rng.nextFloat() <= 0.5) ? perturbLarge(&rng, pathPtr) : perturbSmall(current, &rng, pathPtr);

			const double acceptanceProb = (current.mScalarContrib == 0.0) ? 1.0 : std::min(1.0, proposed.mScalarContrib / current.mScalarContrib);
			if (rng.nextFloat() <= acceptanceProb)
			{
				current = proposed;
			}

			if (current.mScalarContrib != 0.0)
			{
				mCamera2->mFilm->atomicAddSample(current.mRaster, current.mContrib * b / current.mScalarContrib / static_cast<double>(mNumMutations));
			}

			if (mDoVisualize)
			{
				visualizer.drawPath(*pathPtr);
				pathPtr->clear();
			}
		}

		if (mDoVisualize) Viewer::End(visualizer.mFilm);
		FimageIo::Save(visualizer.mFilm->getImage(), "visualized.pfm");
	}

	bool mDoVisualize = true;
	Uint mNumVertices;
	Uint mNumMutations;
	SimpleSubpathSampler2 mSubpathSampler2;
};
