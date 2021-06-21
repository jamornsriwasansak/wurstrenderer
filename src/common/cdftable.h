#pragma once

#include "common/wurst.h"
#include "common/util/file.h"

struct CdfTable
{
	CdfTable() {}

	CdfTable(const std::vector<double> & weights)
	{
		if (weights.size() == 0) return;

		mPdfs.resize(weights.size());
		mCdfs.resize(weights.size());

		double sum = 0.0;
		for (Uint i = 0; i < weights.size(); i++) { sum += std::abs(weights[i]); }

		// normalize weight to get pdf
		for (Uint i = 0; i < weights.size(); i++) { mPdfs[i] = std::abs(weights[i]) / sum; }

		// compute cdf
		mCdfs[0] = mPdfs[0];
		mCdfs[weights.size() - 1] = 1.0;
		for (Uint i = 1; i < weights.size() - 1; i++) { mCdfs[i] = mCdfs[i - 1] + mPdfs[i]; }

		assert(weights.size() != 0);
	}

	bool isEmpty() const
	{
		return mPdfs.empty();
	}

	int sample(double * pdf, double * remappedSample, double sample) const
	{
		auto itor = std::lower_bound(mCdfs.begin(), mCdfs.end(), sample);
		int index = static_cast<int>(std::distance(mCdfs.begin(), itor));
		if (pdf) *pdf = mPdfs[index];
		if (remappedSample) *remappedSample = std::min((index == 0) ? sample / mPdfs[index] : (sample - mCdfs[index - 1]) / mPdfs[index], 1.0 - Math::SmallValue);
		assert(index >= 0 && index < static_cast<int>(mPdfs.size()));
		return index;
	}

	std::vector<double> mPdfs;
	std::vector<double> mCdfs;
};
