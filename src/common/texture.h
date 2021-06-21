#pragma once

#include "common/wurst.h"

template <typename Type>
struct Texture
{
	virtual Type eval(const Vec2 & uv) const = 0;

	// average the whole texture needed for weighting lightsource pdf in cdf table
	virtual Type average() const = 0;
};
