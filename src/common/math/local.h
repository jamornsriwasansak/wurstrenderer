#pragma once

#include "vecmath.h"

struct Local2
{
	static double Cos(const Vec2 & v)
	{
		assert(Math::IsNormalized(v));
		return v[1];
	}

	static double Cos2(const Vec2 & v)
	{
		assert(Math::IsNormalized(v));
		return v[1] * v[1];
	}

	static bool IsFrontface(const Vec2 & v)
	{
		assert(Math::IsNormalized(v));
		return v[1] > 0.0;
	}

	static Vec2 Reflect(const Vec2 & v)
	{
		assert(Math::IsNormalized(v));
		return Vec2(v[0], -v[1]);
	}
};

struct Local3
{
	static double Cos(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return v[1];
	}

	static double Cos2(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return v[1] * v[1];
	}

	static double CosAboutY(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return v[0] / Sin(v);
	}

	static double Cos2AboutY(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return v[0] * v[0] / Sin2(v);
	}

	static double Sin(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return std::sqrt(Sin2(v));
	}

	static double Sin2(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return Math::Clamp(1.0 - Cos2(v), 0.0, 1.0);
	}

	static double SinAboutY(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return v[2] / Sin(v);
	}

	static double Sin2AboutY(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return v[2] * v[2] / Sin2(v);
	}

	static bool IsFrontface(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return v[1] > 0.0;
	}

	static double Tan(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return std::sqrt(Tan2(v));
	}

	static double Tan2(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		double cos2 = Cos2(v);
		double sin2 = Math::Clamp(1.0 - cos2, 0.0, 1.0);
		return sin2 / cos2;
	}

	static Vec3 Reflect(const Vec3 & v)
	{
		assert(Math::IsNormalized(v));
		return Vec3(v[0], -v[1], v[2]);
	}
};