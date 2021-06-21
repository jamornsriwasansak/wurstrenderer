#pragma once

#include <cassert>
#include "common/mapping.h"
#include "scalarmath.h"
#include "vecmath.h"

namespace Math
{
	struct SphericalHarmonics
	{
		static int ShIndexFromLm(const int l, const int m)
		{
			assert(l >= abs(m));
			return l * (l + 1) + m;
		}

		static Ivec2 LmFromShIndex(const int index)
		{
			const int l = int(std::sqrt(index));
			const int m = index - l * (l + 1);
			return Ivec2(l, m);
		}

		static double LegendrePolynomial(const int l, const int m, const double x)
		{
			assert(l >= std::abs(m));
			assert(l <= 25);
			// evaluate an Associated Legendre Polynomial P(l,m,x) at x
			double pmm = 1.0;
			if (m > 0) {
				double somx2 = std::sqrt((1.0 - x)*(1.0 + x));
				double fact = 1.0;
				for (int i = 1; i <= m; i++) {
					pmm *= (-fact) * somx2;
					fact += 2.0;
				}
			}
			if (l == m) return pmm;
			double pmmp1 = x * (2.0*m + 1.0) * pmm;
			if (l == m + 1) return pmmp1;
			double pll = 0.0;
			for (int ll = m + 2; ll <= l; ++ll) {
				pll = ((2.0*ll - 1.0)*x*pmmp1 - (ll + m - 1.0)*pmm) / (ll - m);
				pmm = pmmp1;
				pmmp1 = pll;
			}
			return pll;
		}

		// compute (a - b)! / (a + b)!
		static double DivFactorialFp(const int a, const int b)
		{
			assert(a >= 0);
			assert(b >= 0);
			static double memois[400][400] = { 0.0 };
			if (b == 0) { return 1.0; }
			double fa = static_cast<double>(a);
			double fb = static_cast<double>(b);
			double reciprocal = 1.0;
			for (double x = fa - fb + 1.0; x <= fa + fb; x += 1.0) reciprocal *= x;
			return 1.0 / reciprocal;
		}

		static double LegendreNormalizationK(const int l, const int m)
		{
			// renormalisation constant for SH function
			//double square = ((2.0*l + 1.0)*FactorialFp(l - m)) / (4.0 * Math::Pi * FactorialFp(l + m));
			double square = (2.0*l + 1.0) / (4.0 * Math::Pi) * DivFactorialFp(l, m);
			return std::sqrt(square);
		}

		static double Eval(const int l, const int m, const double theta, const double phi)
		{
			assert(l >= std::abs(m));
			// return a point sample of a Spherical Harmonic basis function
			// l is the band, range [0..N]
			// m in the range [-l..l]
			// phi in the range [0..Pi]
			// theta in the range [0..2*Pi]
			const double sqrt2 = std::sqrt(2.0);
			if (m == 0)
				return LegendreNormalizationK(l, 0) * LegendrePolynomial(l, m, cos(phi));
			else if (m > 0)
				return sqrt2 * LegendreNormalizationK(l, m) * cos(m * theta) * LegendrePolynomial(l, m, cos(phi));
			else // m < 0
				return sqrt2 * LegendreNormalizationK(l, -m) * sin(-m * theta) * LegendrePolynomial(l, -m, cos(phi));
		}

		static double Eval(const int l, const int m, const Vec3 & direction)
		{
			Vec2 thetaPhi = Mapping::SphericalFromWorld(direction);
			return Eval(l, m, thetaPhi[0], thetaPhi[1]);
		}

		static double Eval(const int index, const double theta, const double phi)
		{
			Ivec2 lm = LmFromShIndex(index);
			return Eval(lm[0], lm[1], theta, phi);
		}

		static double Eval(const int index, const Vec3 & direction)
		{
			Ivec2 lm = LmFromShIndex(index);
			return Eval(lm[0], lm[1], direction);
		}
	};
};
