#pragma once

#include "common/sphericalrepresentation.h"

// Spherical Harmonics
template <typename Type>
struct SphHarRep : public SphRep<Type>
{
	static Type DoubleProduct(const SphHarRep & a, const SphHarRep & b)
	{
		assert(a.mCoeffs.size() == b.mCoeffs.size());
		Type result(0.0);
		for (Uint iCoeff = 0; iCoeff < a.mCoeffs.size(); iCoeff++)
		{
			result += a.mCoeffs[iCoeff] * b.mCoeffs[iCoeff];
		}
		return result;
	}

	SphHarRep()
	{}

	SphHarRep(const Uint numCoeffs):
		mCoeffs(numCoeffs)
	{}

	SphHarRep operator*(const double s) const
	{
		SphHarRep result(mCoeffs.size());
		for (Uint iTerm = 0; iTerm < mCoeffs.size(); iTerm++)
		{
			result.mCoeffs[iTerm] = mCoeffs[iTerm] * s;
		}
		return result;
	}

	SphHarRep operator+(const SphHarRep & shr) const
	{
		SphHarRep result(mCoeffs.size());
		for (Uint iCoeff = 0; iCoeff < mCoeffs.size(); iCoeff++)
		{
			result.mCoeffs[iCoeff] = mCoeffs[iCoeff] + shr.mCoeffs[iCoeff];
		}
		return result;
	}

	Type eval(const Vec3 & dir) const override
	{
		Type result;
		for (Uint iCoeff = 0; iCoeff < mCoeffs.size(); iCoeff++)
		{
			Ivec2 lm = Math::SphericalHarmonics::LmFromShIndex((int)iCoeff);
			if (!Math::IsZero(mCoeffs[iCoeff]))
				result += mCoeffs[iCoeff] * Math::SphericalHarmonics::Eval(lm[0], lm[1], dir);
		}
		return result;
	}

	std::vector<Type> mCoeffs;
};