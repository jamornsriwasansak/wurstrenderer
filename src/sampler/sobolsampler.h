#pragma once

#include "common/sampler.h"
#include "ext/sobol/sobol.h"

struct SobolSampler : public Sampler
{
	SobolSampler(const int scramble) : mScramble(scramble)
	{
		throw std::runtime_error("todo: fix");
		mDimension = 0;
		mIndex = 0;
	}

	double get1d() override
	{
		assert(mDimension + 1 < 1024);
		return sobol::sample(mIndex, mDimension++, mScramble);
	}

	Vec2 get2d() override
	{
		assert(mDimension + 2 < 1024);
		Vec2 result;
		result[0] = sobol::sample(mIndex, mDimension++, mScramble);
		result[1] = sobol::sample(mIndex, mDimension++, mScramble);
		return result;
	}

	void nextSample() override
	{
		mDimension = 0;
		mIndex++;
	}

	void setPixel(const Ivec2 & pixel) override
	{
	}

	shared_ptr<Sampler> clone(const int seed) const
	{
		return make_shared<SobolSampler>(mScramble);
	}

	int mIndex;
	int mDimension;
	int mScramble;
};
