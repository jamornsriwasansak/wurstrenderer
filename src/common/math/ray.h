#pragma once

#include <limits>

template <typename Vec>
struct Ray
{
	constexpr static double Epsilon = 1e-4;

	Ray() {}

	Ray(const Vec & o, const Vec & d) : origin(o), direction(d)
	{}

	Ray(const Vec & o, const Vec & d, const double tMin, const double tMax) :
		origin(o), direction(d), tMin(tMin), tMax(tMax)
	{}

	Vec3 p(const double t) const
	{
		return origin + direction * t;
	}

	Vec3 pMax() const
	{
		return origin + direction * tMax;
	}

	Vec origin;
	Vec direction;
	double tMin = Epsilon;
	double tMax = std::numeric_limits<double>::max();
};

using Ray2 = Ray<Vec2>;
using Ray3 = Ray<Vec3>;