#pragma once

#include "vec.h"

using Uint = uint64_t;

using Vec2 = Vec<double, 2, 0>;
using Vec3 = Vec<double, 3, 1>;
using Vec4 = Vec<double, 4, 2>;
using Uvec2 = Vec<Uint, 2, 3>;
using Uvec3 = Vec<Uint, 3, 4>;
using Uvec4 = Vec<Uint, 4, 5>;

/*
inline double ComputeTriangleArea(const Vec3 & a, const Vec3 & b, const Vec3 & c)
{
	Vec3 ab = b - a;
	Vec3 ac = c - a;
	Vec3 abcross = glm::cross(ab, ac);
	return glm::length(abcross) / 2.0f;
}

inline size_t MaxExtent(const Vec3 & v)
{
	if (v.x >= v.y)
	{
		if (v.x >= v.z) { return 0; }
		else { return 2; }
	}
	else // (y > x)
	{
		if (v.y >= v.z) { return 1; }
		else { return 0; }
	}
	assert(false && "code shouldn't have reached this point");
}
*/

namespace Math
{
	// invpi was computed with https://www.wolframalpha.com/input/?i=1%2Fpi
	// according to https://www.jpl.nasa.gov/edu/news/2016/3/16/how-many-decimals-of-pi-do-we-really-need/
	// we should have enough digits to land on the moon.
	const double Pi = 3.141592653589793238462643383279502884197169399375105820974;
	const double InvPi = 0.318309886183790671537767526745028724068919291480912897495;

	template<typename Scalar>
	Scalar Clamp(const Scalar a, const Scalar lower, const Scalar upper)
	{
		return std::max(lower, std::min(upper, a));
	}

	template<typename Scalar>
	Scalar Square(const Scalar a)
	{
		return a * a;
	}

	inline int PositiveMod(const int32_t a, const int32_t b)
	{
		int r = a % b;
		return (r < 0) ? r + b : r;
	}

	template <typename Scalar>
	int FloorToInt(Scalar p)
	{
		return (int)std::floor(p);
	}

	template <typename Scalar> inline int CeilToInt(Scalar p)
	{
		return (int)std::ceil(p);
	}

	inline uint32_t ComputeAlign(uint32_t size, uint32_t align)
	{
		return ((size + align - 1) / align) * align;
	}

	inline uint32_t ComputeAlignPow2(uint32_t size)
	{
		uint32_t value = 1;
		while (value < size) { value *= 2; }
		return value;
	}
}