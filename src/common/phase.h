#pragma once

#include "common/wurst.h"

struct PhaseFunction
{
	virtual double evalPdfW(const MediumVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const = 0;

	virtual Spectrum evalPhase(const MediumVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const = 0;

	virtual Spectrum samplePhaseByPdf(Vec3 * outgoing, double * pdfW, const MediumVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const = 0;
};
