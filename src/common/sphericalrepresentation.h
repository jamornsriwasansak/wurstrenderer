#pragma once

#include "common/wurst.h"

template <typename Type>
struct SphRep
{
	virtual Type eval(const Vec3 & dir) const = 0;

	virtual double pdfA(const Vec3 & dir) const
	{
		throw std::runtime_error("unimpl");
		return 0.0;
	}

	virtual Type sample(Vec3 * dir, double * pdfW, const Vec2 & sample) const
	{
		throw std::runtime_error("unimpl");
		return Type(0.0);
	}
};
