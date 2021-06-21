#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <cmath>

#ifndef __forceinline
#define __forceinline inline
#endif

using Uint = size_t;

// :( fix ambiguous
// double / float
template <typename Scalar>
inline bool isfinite2(Scalar) { return true; }

template <>
inline bool isfinite2(double x) { return std::isfinite(x); }

template <>
inline bool isfinite2(float x) { return std::isfinite(x); }

namespace Math
{
	// invpi was computed with https://www.wolframalpha.com/input/?i=1%2Fpi
	// according to https://www.jpl.nasa.gov/edu/news/2016/3/16/how-many-decimals-of-pi-do-we-really-need/
	// we should have enough digits to land on the moon.
	const double Pi = 3.141592653589793238462643383279502884197169399375105820974;
	const double InvPi = 0.318309886183790671537767526745028724068919291480912897495;
	const double SmallValue = 1e-12;

	template <typename Scalar> inline int CeilToInt(Scalar p)
	{
		return (int)std::ceil(p);
	}

	template <typename Scalar>
	Scalar Clamp(const Scalar a, const Scalar lower, const Scalar upper)
	{
		return std::max(lower, std::min(upper, a));
	}

	uint8_t ClampToUint8(const double a, const uint8_t lower, const uint8_t upper)
	{
		const int ia = int(a);
		const int ilower = int(lower);
		const int iupper = int(upper);
		return uint8_t(std::max(ilower, std::min(ia, iupper)));
	}

	inline Uint ComputeAlign(Uint size, Uint align)
	{
		return ((size + align - 1) / align) * align;
	}

	inline Uint ComputeAlignPow2(Uint size)
	{
		Uint value = 1;
		while (value < size) { value *= 2; }
		return value;
	}

	inline int DivAndCeil(const int a, const int b)
	{
		assert(b != 0);
		return (a + b - 1) / b;
	}

	inline Uint DivAndCeil(const Uint a, const Uint b)
	{
		assert(b != 0);
		return (a + b - 1) / b;
	}

	template <typename Scalar>
	inline Scalar Exp(const Scalar a)
	{
		return std::exp(a);
	}

	inline int Factorial(int k)
	{
		assert(k >= 0);
		static int memois[1000] = { 0 };
		if (k == 0) return 1;
		if (memois[k] != 0) { return memois[k]; }
		const int v = k * Factorial(k - 1);
		memois[k] = v;
		return v;
	}

	inline double FactorialFp(const int k)
	{
		assert(k >= 0);
		static double memois[1000] = { 0.0 };
		if (k == 0) return 1.0;
		if (memois[k] != 0) { return memois[k]; }
		const double v = (double)k * FactorialFp(k - 1);
		memois[k] = v;
		return v;
	}

	template <typename Scalar>
	int FloorToInt(const Scalar p)
	{
		return (int)std::floor(p);
	}

	// similar to python isclose
	inline bool IsApprox(const double a, const double b, const double absTol = SmallValue, const double relTol = SmallValue)
	{
		assert(std::isfinite(a));
		assert(std::isfinite(b));
		const double absA = std::abs(a), absB = std::abs(b);
		const double larger = absA < absB ? absB : absA;
		return std::abs(a - b) <= std::max(larger * relTol, absTol);
	}

	inline bool IsSqrtInteger(const int a)
	{
		if (a < 0) return false;
		int s = static_cast<int>(std::sqrt(static_cast<int>(a)));
		return s * s == a;
	}

	inline bool IsZero(const double a)
	{
		return IsApprox(a, 0.0);
	}

	inline double Length(const double v)
	{
		return v;
	}

	inline int PositiveMod(const int a, const int b)
	{
		int r = a % b;
		return (r < 0) ? r + b : r;
	}

	inline bool SolveQuadratic(double * x0, double * x1, const double a, const double b, const double c)
	{
		if (Math::IsApprox(a, 0))
		{
			const double ans = -c / b;
			if (x0) *x0 = ans;
			if (x1) *x1 = ans;
			return true;
		}

		const double discriminant = b * b - 4 * a * c;
		if (discriminant < 0) return false;

		const double sqrtDiscriminant = std::sqrt(discriminant);
		const double q = (b > 0) ? -0.5 * (b + sqrtDiscriminant) : -0.5 * (b - sqrtDiscriminant);
		if (x0) *x0 = q / a;
		if (x1) *x1 = c / q;
		return true;
	}

	// solve ax^3 + bx^2 + cx + d
	// from https://www.easycalculation.com/algebra/cubic-equation.php with a lot of code cleaning
	// passed the test
	inline int SolveRealCubic(double * x0, double * x1, double * x2, const double ua, const double ub, const double uc, const double ud)
	{
		if (Math::IsApprox(ua, 0)) { return SolveQuadratic(x0, x1, ub, uc, ud) ? 2 : 0; }

		const double b = ub / ua;
		const double c = uc / ua;
		const double d = ud / ua;

		const double b2 = b * b;
		const double bBy3 = (b / 3.0);

		const double q = (3.0 * c - b2) / 9.0;
		const double q3 = q * q * q;
		const double r = -d * 0.5 + b * (9.0 * c - 2.0 * b2) / 54.0;
		const double discriminant = q3 +  r * r;

		if (discriminant > 0)   // One root real, two are complex
		{
			const double sqrtDisc = sqrt(discriminant);
			if (x0) *x0 = -bBy3 + std::cbrt(r + sqrtDisc) + std::cbrt(r - sqrtDisc);
			return 1;
		}
		// The remaining options are all real
		else if (discriminant == 0)  // All roots real, at least two are equal.
		{
			const double r13 = std::cbrt(r);
			if (x0) *x0 = -bBy3 + 2.0*r13;
			if (x1) *x1 = -(r13 + bBy3);
			return 2;
		}
		// Only option left is that all roots are real and unequal (to get here, q < 0)
		else
		{
			const double dum1 = std::acos(r / std::sqrt(-q3));
			const double r13 = 2.0 * std::sqrt(-q);
			if (x0) *x0 = -bBy3 + r13 * std::cos(dum1 / 3.0);;
			if (x1) *x1 = -bBy3 + r13 * std::cos((dum1 + 2.0 * Pi) / 3.0);
			if (x2) *x2 = -bBy3 + r13 * std::cos((dum1 + 4.0 * Pi) / 3.0);
			return 3;
		}
	}

	template<typename Scalar>
	Scalar Square(const Scalar a)
	{
		return a * a;
	}
};
