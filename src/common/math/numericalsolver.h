#pragma once

#include "scalarmath.h"

namespace Math
{
	struct NumericalSolver
	{
		// f1(u, v) = a1(u^2) + a2(v^2) + a3(u) + a4(v) + a5(uv) + a6 = 0
		// f2(u, v) = b1(u^2) + b2(v^2) + b3(u) + b4(v) + b5(uv) + b6 = 0
		static bool SolveTwoQuadraticTwoUnknownsNewton(double * px0, double * px1,
													   const double a1, const double a2, const double a3, const double a4, const double a5, const double a6,
													   const double b1, const double b2, const double b3, const double b4, const double b5, const double b6,
													   const int numSolveIters = 5000, const double epsilon = SmallValue)
		{
			auto f0 = [&](const double u, const double v) { return (a1 * u * u) + (a2 * v * v) + (a3 * u) + (a4 * v) + (a5 * u * v) + a6; };
			auto f1 = [&](const double u, const double v) { return (b1 * u * u) + (b2 * v * v) + (b3 * u) + (b4 * v) + (b5 * u * v) + b6; };

			// create jacobian matrix
			// 
			// J(u,v) = | j11(u,v) j12(u,v) | = | df1(u,v)/du   df1(u,v)/dv |
			//	        | j21(u,v) j22(u,v) |   | df2(u,v)/du   df2(u,v)/dv |

			auto j11 = [&](const double u, const double v) { return 2.0 * a1 * u + a3 + a5 * v; };
			auto j12 = [&](const double u, const double v) { return 2.0 * a2 * v + a4 + a5 * u; };
			auto j21 = [&](const double u, const double v) { return 2.0 * b1 * u + b3 + b5 * v; };
			auto j22 = [&](const double u, const double v) { return 2.0 * b2 * v + b4 + b5 * u; };

			// initial guess
			double x0 = px0 ? *px0 : 1;
			double x1 = px1 ? *px1 : 1;

			for (int iSolve = 0; iSolve < numSolveIters; iSolve++)
			{
				// compute J(x0, x1) = | a b |
				//                     | c d |
				const double a = j11(x0, x1);
				const double b = j12(x0, x1);
				const double c = j21(x0, x1);
				const double d = j22(x0, x1);

				// compute det and early terminate if det is bad
				const double det = a * d - b * c;
				if (Math::IsApprox(det, 0.0)) return false;
				const double invDet = 1.0 / det;

				// compute J^-1(x0, x1)
				const double fx0 = f0(x0, x1);
				const double fx1 = f1(x0, x1);
				const double newX0 = x0 - (invDet) * (fx0 * d - fx1 * b);
				const double newX1 = x1 - (invDet) * (-fx0 * c + fx1 * a);

				if (Math::IsApprox(x0, newX0, epsilon) && Math::IsApprox(x1, newX1, epsilon))
				{
					if (px0) *px0 = x0;
					if (px1) *px1 = x1;
					return true;
				}

				x0 = newX0;
				x1 = newX1;
			}
			return false;
		}
	};
};
