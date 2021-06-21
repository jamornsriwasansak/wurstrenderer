#pragma once

#include "common/wurst.h"

#include <pbrtParser/Scene.h>

struct UtilPbrt
{
	inline static Ivec2 Ivec2FromPbrtVec2i(const pbrt::vec2i & vec2i)
	{
		return Ivec2(vec2i.x, vec2i.y);
	}

	inline static Vec3 Vec3FromPbrtVec3f(const pbrt::vec3f & vec3)
	{
		return Vec3(vec3.x, vec3.y, vec3.z);
	}

	inline static RgbSpectrum RgbSpectrumFromPbrtVec3f(const pbrt::vec3f & vec3)
	{
		return RgbSpectrum(vec3.x, vec3.y, vec3.z);
	}
};
