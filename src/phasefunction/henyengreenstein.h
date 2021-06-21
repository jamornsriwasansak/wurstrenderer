#pragma once

#include "common/wurst.h"

#include "common/path.h"
#include "common/phase.h"

struct HenyeyGreensteinPhase : public PhaseFunction
{
	static inline double PhaseHg(const double cosTheta, const double g)
	{
		double g2 = g * g;
		double denom = 1.0 + g2 + 2 * g * cosTheta;
		return Math::InvPi * 0.25 * (1.0 - g2) / (denom * std::sqrt(denom));
	}

	HenyeyGreensteinPhase(const Spectrum & sigmaS, const double g): mSigmaS(sigmaS), mG(g)
	{
	}

	double evalPdfW(const MediumVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assert(vertex.mPhaseFunction == this);
		throw std::runtime_error("unimplemented");
		return 1.0;
	}

	Spectrum evalPhase(const MediumVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		return Spectrum(PhaseHg(Math::Dot(incoming, outgoing), mG));
	}

	Spectrum samplePhaseByPdf(Vec3 * outgoing, double * pdfW, const MediumVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const override
	{
		// this is sampling code from pbrt
		// haven't tested yet
		double cosTheta;
		if (std::abs(mG) < 1e-3)
			cosTheta = 1.0 - 2.0 * sample[0];
		else {
			double sqrTerm = (1.0 - mG * mG) / (1.0 - mG + 2.0 * mG * sample[0]);
			cosTheta = (1.0 + mG * mG - sqrTerm * sqrTerm) / (2.0 * mG);
		}

		double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
		double phi = 2.0 * Math::Pi * sample[1];
		Vec3 v1, v2;
		CoordFrame3 basis(incoming);
		*outgoing = basis.toWorld(Mapping::SphereFromSquare(Vec2(sinTheta, cosTheta)));
		return PhaseHg(-cosTheta, mG) * mSigmaS;
	}

	double mG;
	Spectrum mSigmaS;
};
