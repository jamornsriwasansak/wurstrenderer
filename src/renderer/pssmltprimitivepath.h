#pragma once


#include "common/wurst.h"

#include "common/camera.h"
#include "common/mlt.h"
#include "common/parallel.h"
#include "common/renderer.h"
#include "common/scene.h"
#include "common/sampler.h"

#include "renderer/primitivepath.h"

struct PssmltPrimitivePathRenderer : public SingleCameraRenderer
{
	PssmltPrimitivePathRenderer(shared_ptr<const Camera> & camera,
								shared_ptr<const Scene> & scene,
								shared_ptr<Sampler> & sampler,
								const int mNumPathVertices,
								const int numSamples) :
		SingleCameraRenderer(camera, scene, sampler),
		mNumVertices(mNumPathVertices),
		mNumSamples(numSamples * Math::Volume(camera->mFilm->mResolution)),
		mBaseRenderer(camera, scene, sampler, mNumPathVertices, 0, false)
	{
	}

	// note: changed on 6 July 2019. Untested
	PssmltState perturbAndEval(const PssmltState & prev, const bool isLargeStep, Rng * rng) const
	{
		PssmltState result = (isLargeStep) ? PssmltState::PerturbPrimarySampleLarge(prev, rng) : PssmltState::PerturbPrimarySampleSmall(prev, rng);
		PssmltSampler pssmltSampler(&result.mPrimarySample);
		Splat splat;
		mBaseRenderer.evalContribution(&splat, &pssmltSampler);
		result.mContrib = splat.mContrib;
		result.mRaster = splat.mRaster;
		result.mScalarContrib = Math::Length(result.mContrib);
		return result;
	}

	void render() override
	{
		// scale the whole image down by numsamples
		mCamera->mFilm->mAtomicScalingFactor.store(1.0 / static_cast<double>(mNumSamples));

		Rng rng(0);

		int numNormalizationFactorSamples = 0;

		// estimate normalization factor
		double b = 0.0;
		PssmltState state(mNumVertices * 2);
		for (int iB = 0; iB < 1000000; iB++)
		{
			b += perturbAndEval(state, true, &rng).mScalarContrib;
			numNormalizationFactorSamples += 1;
		}

		int numDesiredMutationPerWorkUnit = 200000;

		std::mutex mutex;
		int numCurrentSamples = 0;

		Parallel::Split(0, mNumSamples, numDesiredMutationPerWorkUnit, [&](const int start, const int end)
		{
			Rng rng(start);
			PssmltState current = perturbAndEval(state, true, &rng);
			int localNumLargeStep = 0;
			double localB = 0.0;
			for (int iMutation = start; iMutation < end; iMutation++)
			{
				const bool doPerturbLarge = (rng.nextFloat() < 0.3);
				PssmltState proposed = perturbAndEval(current, doPerturbLarge, &rng);
				const double acceptanceProb = (current.mScalarContrib == 0.0) ? 1.0 : std::min(1.0, proposed.mScalarContrib / current.mScalarContrib);

				if (doPerturbLarge)
				{
					localNumLargeStep += 1;
					localB += proposed.mScalarContrib;
				}

				// splat proposed
				if (proposed.mScalarContrib != 0.0)
				{
					mCamera->mFilm->atomicAddSample(proposed.mRaster, proposed.mContrib * acceptanceProb / proposed.mScalarContrib);
				}

				// splat current
				if (current.mScalarContrib != 0.0)
				{
					mCamera->mFilm->atomicAddSample(current.mRaster, current.mContrib * (1.0 - acceptanceProb) / current.mScalarContrib);
				}

				if (rng.nextFloat() <= acceptanceProb)
				{
					current = proposed;
				}

				// update film result
				if ((iMutation + 1) % 5000000 == 0)
				{
					mutex.lock();
					numCurrentSamples += 5000000;
					b += localB;
					numNormalizationFactorSamples += localNumLargeStep;
					mCamera->mFilm->mAtomicScalingFactor.store(b / static_cast<double>(numNormalizationFactorSamples * numCurrentSamples));

					if (numCurrentSamples % 100000000 == 0)
					{
						std::cout << "write" << std::endl;
						std::string filename = "result_" + std::to_string(numCurrentSamples) + ".pfm";
						FimageIo::Save(mCamera->mFilm->getImage(), filename);
					}

					mutex.unlock();
					Viewer::Redraw(mCamera->mFilm);
				}
			}
		});
	}

	int mNumVertices;
	int mNumSamples;
	PrimitivePathRenderer mBaseRenderer;
};
