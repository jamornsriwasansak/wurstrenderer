#pragma once

#include "common/wurst.h"
#include "common/rng.h"

struct Sampler
{
	Sampler() {}

	virtual double get1d() = 0;

	virtual Vec2 get2d() = 0;

	virtual void nextSample() { throw std::runtime_error("unimpl"); }

	virtual void setPixel(const Ivec2 & pixel) { throw std::runtime_error("unimpl"); }

	virtual shared_ptr<Sampler> clone(const int seed) const { throw std::runtime_error("unimpl"); return nullptr; }
};
