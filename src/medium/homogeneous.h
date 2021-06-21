#pragma once

#include "common/wurst.h"

#include "common/medium.h"
#include "common/path.h"
#include "phasefunction/isotropic.h"

struct HomogeniusMedium : public Medium
{
	static shared_ptr<const HomogeniusMedium> FromSimple2Json(const json & json)
	{
		const RgbSpectrum sigmaA = Util::Json::RgbFromJson(json, "sigma_a");
		const RgbSpectrum sigmaS = Util::Json::RgbFromJson(json, "sigma_s");

		shared_ptr<const PhaseFunction> phase = make_shared<const IsotropicPhase>(sigmaS);

		return make_shared<const HomogeniusMedium>(sigmaS, sigmaA, phase);
	}

	HomogeniusMedium(const Spectrum sigmaS, const Spectrum sigmaA, const shared_ptr<const PhaseFunction> & phaseFunction):
		mSigmaT(sigmaS + sigmaA),
		mPhaseFunction(phaseFunction)
	{
	}

	Spectrum evalTransmittance(const double dist) const override
	{
		return Math::Exp(-mSigmaT * dist);
	}

	shared_ptr<MediumVertex> sampleMediumVertex(Spectrum * s, const Ray3 & ray, const Vec2 & sample) const override
	{
		assert(Math::IsNormalized(ray.direction));

		int channel = std::min(int(sample[0] * Spectrum::NumElements), int(Spectrum::NumElements - 1));
		const double dist = -std::log(1.0 - sample[1]) / mSigmaT[channel];
		const double t = std::min(dist, ray.tMax);
		const bool sampledMedium = t < ray.tMax;

		// compute tr and pdf
		const Spectrum tr = Spectrum(Math::Exp(-mSigmaT * t));
		const Spectrum density = sampledMedium ? mSigmaT * tr : tr;
		double pdf = 0.0;
		for (int i = 0; i < Spectrum::NumElements; i++) { pdf += density[i]; }
		pdf /= double(Spectrum::NumElements);
		pdf = (pdf == 0.0) ? 1.0 : pdf;

		// assign values and vertex information
		if (s) *s = tr / pdf;
		if (sampledMedium)
		{
			shared_ptr<MediumVertex> vertex = make_shared<MediumVertex>();
			vertex->mCoordFrame = CoordFrame3(-ray.direction);
			vertex->mGeometryNormal = -ray.direction;
			vertex->mPhaseFunction = mPhaseFunction.get();
			vertex->mPosition = ray.p(t);
			vertex->mMediumFrontface = this;
			vertex->mMediumBackface = this;
			return vertex;
		}
		else
		{
			return nullptr;
		}
	}

	Spectrum mSigmaT;
	shared_ptr<const PhaseFunction> mPhaseFunction;
};
