#pragma once

#include "vecprim.h"

#if __AVX2__
#ifndef AVX2
#define AVX2
#endif
#endif

#ifdef AVX2
#include "vecsimd.h"
#endif

#ifdef AVX2
using Vec4 = Math::VecSimd4<4, 2>;
using Vec3 = Math::VecSimd4<3, 1>;
#else
using Vec4 = Math::VecPrim<double, 4, 2>;
using Vec3 = Math::VecPrim<double, 3, 1>;
#endif
using Vec2 = Math::VecPrim<double, 2, 0>;

using Uvec4 = Math::VecPrim<Uint, 4, 5>;
using Uvec3 = Math::VecPrim<Uint, 3, 4>;
using Uvec2 = Math::VecPrim<Uint, 2, 3>;

using Ivec4 = Math::VecPrim<int, 4, 8>;
using Ivec3 = Math::VecPrim<int, 3, 7>;
using Ivec2 = Math::VecPrim<int, 2, 6>;

// needed so that we can convert between Uvec, Ivec and SimdVec easily
	
namespace Math
{
#ifdef AVX2
	template<typename Scalar, size_t Size, size_t Id>
	template<size_t OtherSize, size_t OtherId>
	inline VecPrim<Scalar, Size, Id>::VecPrim(const VecSimd4<OtherSize, OtherId> & u)
	{
		static_assert(Size == OtherSize, "Vec: Must initialize with the vector of equal Size");
		for (size_t i = 0; i < Size; i++) { this->v[i] = Scalar(u[i]); }
		assert(IsFinite(*this));
	}

	template<size_t Size, size_t Id>
	template<typename OtherScalar, size_t OtherSize, size_t OtherId>
	inline VecSimd4<Size, Id>::VecSimd4(const VecPrim<OtherScalar, OtherSize, OtherId>& ov)
	{
		static_assert(Size == OtherSize, "VecSimd4: Size mismatch");
		for (size_t i = 0; i < Size; i++) { u[i] = double(ov[i]); }
		assert(IsFinite(*this));
	}
#endif

	template<typename VecType>
	inline VecType BarycentricInterpolate(const VecType & a, const VecType & b, const VecType & c, const Vec2 & barycentric)
	{
		return (VecType(1) - barycentric[0] - barycentric[1]) * a + barycentric[0] * b + barycentric[1] * c;
	}

	template<typename VecType>
	inline VecType FlipY(const VecType & a)
	{
		VecType result = a;
		result[1] = -result[1];
		return result;
	}

	template<typename VecType>
	inline VecType Interpolate(const VecType & a, const VecType & b, const double t)
	{
		return a * t + (1 - t) * b;
	}

	template<typename VecType>
	inline bool IsFrontface(const VecType & direction, const VecType & normal)
	{
		return Math::Dot(direction, normal) > 0;
	}

	inline double ParallelogramArea(const Vec3 & a, const Vec3 & b)
	{
		// area = |a x b|
		// area = |a||b|sin(theta)
		// area^2 = |a|^2 |b|^2 sin^2(theta)
		// area^2 = |a|^2 |b|^2 (1.0 - cos^2(theta))

		double area2 = Math::Length2(a) * Math::Length2(b) * (1.0 - Math::Square(Math::Dot(a, b)));
		return std::sqrt(area2);
	}

	template<typename Vector>
	Vector RadiansFromDegrees(const Vector & a)
	{
		return a * Math::Pi / 180.0;
	}

	inline Vec3 Reflect(const Vec3 & i, const Vec3 & n)
	{
		return i - 2.0 * Math::Dot(n, i) * n;
	}
}