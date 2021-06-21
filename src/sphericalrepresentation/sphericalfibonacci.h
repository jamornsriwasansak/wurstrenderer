#pragma once

#include "common/sphericalrepresentation.h"

template <typename Type>
struct SphFiboRep : public SphRep<Type>
{
	static Spectrum DoubleProduct(const SphFiboRep & a, const SphFiboRep & b)
	{
		assert(a.mData.size() == b.mData.size());
		Spectrum result(0.0);
		for (Uint iCoeff = 0; iCoeff < a.mData.size(); iCoeff++)
		{
			result += a.mData[iCoeff] * b.mData[iCoeff];
		}
		return result;
	}

	SphFiboRep()
	{}

	SphFiboRep(const Uint numCoeffs):
		mData(numCoeffs, Type(0.0))
	{}

	SphFiboRep(const std::vector<Type> & coeffs):
		mData(coeffs)
	{}

	Type eval(const Vec3 & direction) const override
	{
		std::array<FiboInterpolateInfo, 9> interp; 
		Mapping::SphFiboIndicesFromWorld(&interp, direction, static_cast<int>(mData.size()));
		Type result(0.0);
		for (int i = 0; i < 9; i++)
		{
			result += mData[interp[i].index] * interp[i].weight;
		}
		return result;
	}

	std::vector<Type> mData;
};
