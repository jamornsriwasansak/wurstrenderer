#pragma once

#include "common/sampler.h"

struct RandomSampler : public Sampler
{
	RandomSampler(const int seed):
		mRng(seed)
	{
	}

	double get1d() override
	{
		return mRng.nextFloat();
	}

	Vec2 get2d() override
	{
		return Vec2(mRng.nextFloat(), mRng.nextFloat());
	}

	void nextSample() override
	{
	}

	void setPixel(const Ivec2 & pixel) override
	{
	}

	shared_ptr<Sampler> clone(const int seed) const override
	{
		shared_ptr<Sampler> result = make_shared<RandomSampler>(seed);
		return result;
	}

	Rng mRng;
};
