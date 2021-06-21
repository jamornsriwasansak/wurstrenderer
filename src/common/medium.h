#pragma once

#include "common/wurst.h"

#include "common/path.h"
#include "common/phase.h"

struct Medium
{
	virtual Spectrum evalTransmittance(const double dist) const = 0;

	virtual shared_ptr<MediumVertex> sampleMediumVertex(Spectrum * s, const Ray3 & ray, const Vec2 & sample) const = 0;
};
