#pragma once

#include "common/bsdf.h"
#include "common/path.h"
#include "common/texture.h"

// refer to pbrt v3. page 534
struct OrenNayarBsdf : public Bsdf
{
	OrenNayarBsdf(const shared_ptr<const Texture<Spectrum>> diffuseReflectance, const shared_ptr<const Texture<Spectrum>> stdDeviation):
		mDiffuseReflectance(diffuseReflectance), mStdDeviation(stdDeviation)
	{}

	Spectrum evalBsdf(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));

		if ((Local3::Cos(incoming) <= 0.0) || (Local3::Cos(outgoing) <= 0.0)) { return Spectrum(0.0); }

		double sinPhiI = Local3::Sin(incoming);
		double sinPhiO = Local3::Sin(outgoing);

		double maxCos = 0.0;
		double sinAlpha = 0.0;
		double tanBeta = 0.0;
		if ((Local3::Sin(incoming) > 1e-4) && (Local3::Sin(outgoing) > 1e-4))
		{
			double sinThetaI = Local3::SinAboutY(incoming);
			double cosThetaI = Local3::CosAboutY(incoming);
			double sinThetaO = Local3::SinAboutY(outgoing);
			double cosThetaO = Local3::CosAboutY(outgoing);
			double dCos = cosThetaI * cosThetaO + sinThetaI * sinThetaO;
			maxCos = std::max(0.0, dCos);

			if (std::abs(Local3::Cos(incoming)) > std::abs(Local3::Cos(outgoing)))
			{
				sinAlpha = sinPhiO;
				tanBeta = sinPhiI / Local3::Cos(incoming);
			}
			else
			{
				sinAlpha = sinPhiI;
				tanBeta = sinPhiO / Local3::Cos(outgoing);
			}
		}

		const double stdDeviation = mStdDeviation->eval(vertex.mTextureCoord)[0];
		assert(stdDeviation > Math::SmallValue);
		const double stdDeviation2 = stdDeviation * stdDeviation;
		const double A = 1.0 - stdDeviation2 / (2 * (stdDeviation2 + 0.33));
		const double B = 0.45 * stdDeviation2 / (stdDeviation2 + 0.09);

		return mDiffuseReflectance->eval(vertex.mTextureCoord) * Math::InvPi * (A + B * maxCos * sinAlpha * tanBeta);
	}

	double evalPdfW(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));
		if (Local3::Cos(outgoing) <= 0.0) { return 0.0; }
		if (Local3::Cos(incoming) <= 0.0) { return 0.0; }
		return Local3::Cos(outgoing) * Math::InvPi;
	}

	BsdfType getType() const override
	{
		return BsdfType::Reflect;
	}

	Spectrum sampleBsdfCosByPdf(Vec3 * outgoing, double * pdfW, const SurfaceVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const override
	{
		assert(Math::IsNormalized(incoming));
		if (Local3::Cos(incoming) < 0.0)
		{
			if (pdfW) *pdfW = 0.0;
			return Spectrum(0.0);
		}
		Vec3 out = Mapping::CosineWeightedHemisphereFromSquare(sample);
		if (outgoing) *outgoing = out;
		if (pdfW) *pdfW = Local3::Cos(out) * Math::InvPi;
		return evalBsdf(vertex, out, incoming) * Math::Pi;
	}

	Spectrum totalEnergy(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		if (!Local3::IsFrontface(incoming) || !Local3::IsFrontface(outgoing)) return Spectrum(0.0);
		return mDiffuseReflectance->eval(vertex.mTextureCoord);
	}

	shared_ptr<const Texture<Spectrum>> mDiffuseReflectance;
	shared_ptr<const Texture<Spectrum>> mStdDeviation;
};

