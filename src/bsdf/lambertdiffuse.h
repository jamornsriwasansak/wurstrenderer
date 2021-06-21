#pragma once

#include "common/bsdf.h"
#include "common/path.h"
#include "common/texture.h"

struct LambertDiffuseBsdf : public Bsdf
{
	LambertDiffuseBsdf(shared_ptr<const Texture<Spectrum>> diffuseReflectance):
		mDiffuseReflectance(diffuseReflectance)
	{
	}

	// 3d stuffs

	Spectrum evalBsdf(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));
		if (Local3::Cos(incoming) <= 0.0) { return Spectrum(0.0); }
		if (Local3::Cos(outgoing) <= 0.0) { return Spectrum(0.0); }

		Spectrum diffRefl = mDiffuseReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(diffRefl, Spectrum(1.0)));

		return diffRefl * Math::InvPi;
	}

	double evalPdfW(const SurfaceVertex & , const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));
		if (Local3::Cos(incoming) <= 0.0) { return 0.0; }
		if (Local3::Cos(outgoing) <= 0.0) { return 0.0; }
		return Local3::Cos(outgoing) * Math::InvPi;
	}

	Spectrum sampleBsdfCosByPdf(Vec3 * outgoing, double * pdfW, const SurfaceVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const override
	{
		assert(Math::IsNormalized(incoming));
		if (Local3::Cos(incoming) <= 0.0) { return Spectrum(0.0); }
		Vec3 out = Mapping::CosineWeightedHemisphereFromSquare(sample);
		const double cosine = Local3::Cos(out);
		if (outgoing) *outgoing = out;
		if (pdfW) *pdfW = std::max(cosine, 0.0) * Math::InvPi;
		if (cosine <= 0.0) { return Spectrum(0.0); }

		Spectrum diffRefl = mDiffuseReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(diffRefl, Spectrum(1.0)));
		return diffRefl;
	}

	Spectrum totalEnergy(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		if (!Local3::IsFrontface(incoming) || !Local3::IsFrontface(outgoing)) return Spectrum(0.0);
		return mDiffuseReflectance->eval(vertex.mTextureCoord);
	}
	
	// 2d stuffs

	Spectrum evalBsdf2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));
		if (Local2::Cos(incoming) <= 0.0) { return Spectrum(0.0); }
		if (Local2::Cos(outgoing) <= 0.0) { return Spectrum(0.0); }

		Spectrum diffRefl = mDiffuseReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(diffRefl, Spectrum(1.0)));
		return diffRefl * 0.5;
	}

	double evalPdfT2(const SurfaceVertex & , const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));
		if (Local2::Cos(incoming) <= 0.0) { return 0.0; }
		if (Local2::Cos(outgoing) <= 0.0) { return 0.0; }
		return Local2::Cos(outgoing) * 0.5;
	}

	Spectrum sampleBsdfCosByPdf2(Vec2 * outgoing, double * pdfT, const SurfaceVertex & vertex, const Vec2 & incoming, const double sample) const override
	{
		assert(Math::IsNormalized(incoming));
		if (Local2::Cos(incoming) <= 0.0) { return Spectrum(0.0); }
		Vec2 out = Mapping::CosineWeightedHemicircleFromSegment(sample);
		const double cosine = Local2::Cos(out);
		if (outgoing) *outgoing = out;
		if (pdfT) *pdfT = Math::Pi;
		if (cosine <= 0.0) { return Spectrum(0.0); }

		Spectrum diffRefl = mDiffuseReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(diffRefl, Spectrum(1.0)));
		return diffRefl;
	}

	Spectrum totalEnergy2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		if (!Local2::IsFrontface(incoming) || !Local2::IsFrontface(outgoing)) return Spectrum(0.0);
		return mDiffuseReflectance->eval(vertex.mTextureCoord);
	}

	// common stuffs

	BsdfType getType() const override
	{
		return BsdfType::Reflect;
	}

	shared_ptr<const Texture<Spectrum>> mDiffuseReflectance;
};
