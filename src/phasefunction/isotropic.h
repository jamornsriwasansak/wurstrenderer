#pragma once

#include "common/wurst.h"

#include "common/path.h"
#include "common/phase.h"

struct IsotropicPhase : public PhaseFunction
{
	IsotropicPhase(const Spectrum & sigmaS): mSigmaS(sigmaS)
	{
	}

	double evalPdfW(const MediumVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		return Math::InvPi * 0.25;
	}

	Spectrum evalPhase(const MediumVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		return Math::InvPi * 0.25 * mSigmaS;
	}

	Spectrum samplePhaseByPdf(Vec3 * outgoing, double * pdfW, const MediumVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const override
	{
		assert(vertex.mPhaseFunction == this);
		if (outgoing) *outgoing = Mapping::SphereFromSquare(sample);
		if (pdfW) *pdfW = Math::InvPi * 0.25;
		return mSigmaS;
	}

	Spectrum mSigmaS;
};
