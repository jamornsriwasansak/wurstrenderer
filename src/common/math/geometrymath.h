#pragma once

#include "coordframe.h"
#include "scalarmath.h"
#include "vecmath.h"

namespace Math
{
	inline double ComputeTriangleArea(const Vec3 & a, const Vec3 & b, const Vec3 & c)
	{
		Vec3 ab = b - a;
		Vec3 ac = c - a;
		Vec3 abcross = Math::Cross(ab, ac);
		return Math::Length(abcross) / 2.0f;
	}

	inline bool IntersectSphere(double * t0, double * t1, const Vec3 & center, const double radius, const Vec3 & rayOrigin, const Vec3 & rayDirection)
	{
		assert(t0 != nullptr);
		assert(t1 != nullptr);
		const Vec3 l = rayOrigin - center;
		const double a = Math::Dot(rayDirection, rayDirection);
		const double b = 2.0 * Math::Dot(rayDirection, l);
		const double c = Math::Dot(l, l) - radius * radius;
		if (!SolveQuadratic(t0, t1, a, b, c)) return false;
		if (*t0 > *t1) std::swap(*t0, *t1);
		return true;
	}

	inline bool IntersectPlaneSegment(double * t0, double planeZ, const Vec3 & p1, const Vec3 & p2)
	{
		assert(t0 != nullptr);
		*t0 = (planeZ - p1[2]) / (p2[2] - p1[2]);
		return (*t0 >= 0 && *t0 <= 1);
	}

	/*
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
	}
	*/
}