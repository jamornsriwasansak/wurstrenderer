#pragma once

#include "vecmath.h"

struct CoordFrame2
{
	CoordFrame2() {}

	CoordFrame2(const Vec2 & normal, const bool doCheckNormalize = true)
	{
		if (doCheckNormalize) { assert(Math::IsNormalized(normal)); }
		y = normal;
		x = Vec2(normal[1], -normal[0]);
	}

	CoordFrame2(const Vec2 & xBasis, const Vec2 & yBasis, const bool doCheckNormalize = true):
		x(xBasis),
		y(yBasis)
	{
		if (doCheckNormalize)
		{
			assert(Math::IsNormalized(xBasis));
			assert(Math::IsNormalized(yBasis));
			assert(Math::IsApprox(Math::Dot(xBasis, yBasis), 0.0));
		}
	}

	Vec2 toLocal(const Vec2 & global) const
	{
		return Vec2(Math::Dot(global, x), Math::Dot(global, y));
	}

	Vec2 toWorld(const Vec2 & local) const
	{
		return local[0] * x + local[1] * y;
	}

	Vec2 x, y;
};

struct CoordFrame3
{
	CoordFrame3():
		x(Vec3(1.0, 0.0, 0.0)),
		y(Vec3(0.0, 1.0, 0.0)),
		z(Vec3(0.0, 0.0, 1.0))
	{
	}

	CoordFrame3(const Vec3 & normal, const bool doAssertNormalize = true)
	{
		if (doAssertNormalize) { assert(Math::IsNormalized(normal)); }
		y = normal;
		// Building an Orthonormal Basis, Revisited - Pixar Graphics 2017
		double sign = std::copysign((double)(1.0), normal[1]);
		const double a = -1.0f / (sign + normal[1]);
		const double b = normal[2] * normal[0] * a;
		x = Vec3(sign + normal[0] * normal[0] * a, -normal[0], b);
		z = Vec3(sign * b, -sign * normal[2], 1.0f + sign * normal[2] * normal[2] * a);
	}

	CoordFrame3(const Vec3 & xBasis, const Vec3 & yBasis, const Vec3 & zBasis, const bool doAssertNormalize = true):
		x(xBasis),
		y(yBasis),
		z(zBasis)
	{
		if (doAssertNormalize)
		{
			assert(Math::IsNormalized(xBasis));
			assert(Math::IsNormalized(yBasis));
			assert(Math::IsNormalized(zBasis));
			assert(Math::IsApprox(Math::Dot(Math::Cross(xBasis, yBasis), zBasis), 1.0));

			assert(Math::IsApprox(Math::Dot(xBasis, yBasis), 0.0));
			assert(Math::IsApprox(Math::Dot(yBasis, zBasis), 0.0));
			assert(Math::IsApprox(Math::Dot(xBasis, zBasis), 0.0));
		}
	}

	Vec3 toLocal(const Vec3 & global) const
	{
		return Vec3(Math::Dot(global, x), Math::Dot(global, y), Math::Dot(global, z));
	}

	Vec3 toWorld(const Vec3 & local) const
	{
		return local[0] * x + local[1] * y + local[2] * z;
	}

	Vec3 x, y, z;
};