#pragma once

#include "common/wurst.h"

struct Rng
{
	Rng() : mGenerator64(std::time(0)) {}

	Rng(const Uint seed) : mGenerator64(seed) {}

	double nextFloat(const double start, const double end) 
	{
		assert(end >= start);
		return std::uniform_real_distribution<double>(start, end)(mGenerator64);
	}

	double nextFloat() 
	{
		return this->nextFloat(0, 1);
	}

	int nextInt(const int start, const int end)
	{
		assert(end >= start);
		return std::uniform_int_distribution<int>(start, end)(mGenerator64);
	}

	std::mt19937_64 mGenerator64;
};
