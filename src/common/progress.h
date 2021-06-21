#pragma once

#include "common/wurst.h"

struct ProgressReport
{
	ProgressReport(const int64_t totalNumSamples):
		mStartTime(std::chrono::system_clock::now()),
		mTotalNumSamples(totalNumSamples),
		mNumSamples(0)
	{
	}

	ProgressReport(const int64_t numSamples, const int64_t numPixels): ProgressReport(numSamples * numPixels)
	{
	}

	double getPercent() const
	{
		return double(mNumSamples) / double(mTotalNumSamples);
	}

	int64_t getTimeMilliSecs() const
	{
		return int64_t(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - mStartTime).count());
	}

	int64_t getEtcMilliSecs() const
	{
		int64_t elapsedTime = getTimeMilliSecs();
		int64_t estimateTotalTime = elapsedTime * (mTotalNumSamples + 1) / (mNumSamples + 1);
		return estimateTotalTime - elapsedTime;
	}

	void increment(const int64_t numIncreasedSamples)
	{
		mNumSamples += numIncreasedSamples;
	}

	std::chrono::system_clock::time_point mStartTime;
	std::chrono::system_clock::time_point mLastWrite;
	std::mutex mutex;

	int64_t mTotalNumSamples;
	filesystem::path mPath;
	std::atomic_int64_t mNumSamples;
	shared_ptr<const Film> mFilm;
};
