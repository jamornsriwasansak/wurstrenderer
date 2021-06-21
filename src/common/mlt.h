#pragma once

struct MltUtil
{
	static double PerturbExpLog(const double value, const double s1, const double s2, double r)
	{
		double result;
		if (r < 0.5) {
			r = r * 2.0;
			result = value + s2 * exp(-log(s2 / s1) * r);
			if (result > 1.0)
				result -= 1.0;
		}
		else {
			r = (r - 0.5) * 2.0;
			result = value - s2 * exp(-log(s2 / s1) * r);
			if (result < 0.0)
				result += 1.0;
		}
		return result;
	}
};

// TODO:: implement lazy eval

struct PssmltState
{
	__forceinline static PssmltState PerturbPrimarySampleLarge(const PssmltState & prev, Rng * rng)
	{
		PssmltState result(prev.mDimension);
		for (Uint i = 0; i < prev.mDimension; i++)
		{
			result.mPrimarySample[i] = rng->nextFloat();
		}
		return result;
	}

	__forceinline static PssmltState PerturbPrimarySampleSmall(const PssmltState & prev, Rng * rng)
	{
		PssmltState result(prev.mDimension);
		for (Uint i = 0; i < prev.mDimension; i++)
		{
			result.mPrimarySample[i] = MltUtil::PerturbExpLog(prev.mPrimarySample[i], 1.0 / 1024.0, 1.0 / 64.0, rng->nextFloat());
		}
		return result;
	}

	PssmltState(const Uint dimension): mDimension(dimension), mPrimarySample(dimension)
	{
	}

	Uint mDimension;
	Spectrum mContrib = Spectrum(0.0);
	Vec2 mRaster;
	std::vector<double> mPrimarySample;
	double mScalarContrib;
};

struct PssmltSampler : public Sampler
{
	PssmltSampler(const std::vector<double> * rnds) :
		mCounter(0),
		mRnds(rnds)
	{
	}

	double get1d() override
	{
		mCounter += 1;
		assert(mCounter <= mRnds->size());
		return (*mRnds)[mCounter - 1];
	}

	Vec2 get2d() override
	{
		mCounter += 2;
		assert(mCounter <= mRnds->size());
		return Vec2((*mRnds)[mCounter - 2], (*mRnds)[mCounter - 1]);
	}

	const std::vector<double> * mRnds;
	Uint mCounter;
};
