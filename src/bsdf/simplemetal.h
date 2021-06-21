#pragma once

#include "common/wurst.h"

#include "common/bsdf.h"
#include "common/path.h"
#include "common/texture.h"

struct SimpleMetalBsdf : public Bsdf
{
	static double AngleFromPhongExponent(const double exponent)
	{
		// TODO:: find better fit than 0.5 and 1
		const double base = 0.5;
		const double powerScale = 1;
		return std::pow(base, powerScale * exponent) * Math::Pi;
	}

	SimpleMetalBsdf(const shared_ptr<const Texture<Spectrum>> & specularReflectance, const shared_ptr<const Texture<double>> & specularLobe) :
		mReflectance(specularReflectance),
		mExponent(specularLobe)
	{
	}

	// 3d stuffs

	Spectrum evalBsdf(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));
		if ((Local3::Cos(incoming) <= 0.0) || (Local3::Cos(outgoing) <= 0.0)) { return Spectrum(0.0); }

		const double exponent = mExponent->eval(vertex.mTextureCoord);
		const double halfAngle = AngleFromPhongExponent(exponent);
		const double cosHalfAngle = std::cos(halfAngle);

		const Vec3 reflect = Local3::Reflect(-incoming);
		if (Local3::Cos(reflect) <= 0.0f) { return Spectrum(0.0); }
		const double cosReflect = Math::Dot(reflect, outgoing);
		if (cosReflect < cosHalfAngle) { return Spectrum(0.0); }

		Spectrum specRefl = mReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(specRefl, Spectrum(1.0)));

		return specRefl * 0.5 * Math::InvPi / (1.0 - cosHalfAngle);
	}

	double evalPdfW(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));

		// the logic here is that sample function still have a chance to generate direction under the surface.
		// we still have to evaluate the function if incoming / outgoing goes under the surface.
		if ((Local3::Cos(incoming) <= 0.0) && (Local3::Cos(outgoing) <= 0.0)) { return 0.0; }

		const double exponent = mExponent->eval(vertex.mTextureCoord);
		const double halfAngle = AngleFromPhongExponent(exponent);
		const double cosHalfAngle = std::cos(halfAngle);
		const Vec3 reflect = Local3::Reflect(-incoming);
		const double cosReflect = Math::Dot(reflect, outgoing);
		if (cosReflect < cosHalfAngle) { return 0.0; }

		return 0.5 * Math::InvPi / (1.0 - cosHalfAngle);
	}

	Spectrum sampleBsdfCosByPdf(Vec3 * outgoing, double * pdfW, const SurfaceVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const override
	{
		assert(Math::IsNormalized(incoming));

		if (Local3::Cos(incoming) <= 0.0) { return Spectrum(0.0); }

		const double exponent = mExponent->eval(vertex.mTextureCoord);
		const double halfAngle = AngleFromPhongExponent(exponent);
		const Vec3 reflect = Local3::Reflect(-incoming);

		const CoordFrame3 reflectBasis(reflect);
		const Vec3 out = reflectBasis.toWorld(Mapping::SolidAngleFromSquare(sample, halfAngle));
		const double cosReflect = Math::Dot(out, reflect);
		const double cosHalfAngle = std::cos(halfAngle);

		if (outgoing) *outgoing = out;
		if (pdfW) *pdfW = 0.5 * Math::InvPi / (1.0 - cosHalfAngle);

		if (Local3::Cos(reflect) <= 0.0 || Local3::Cos(out) <= 0.0 || cosReflect <= cosHalfAngle) { return Spectrum(0.0); }
		
		Spectrum specRefl = mReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(specRefl, Spectrum(1.0)));

		return specRefl * Local3::Cos(out);
	}

	Spectrum totalEnergy(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		if (!Local3::IsFrontface(incoming) || !Local3::IsFrontface(outgoing)) return Spectrum(0.0);
		return mReflectance->eval(vertex.mTextureCoord);
	}

	// 2d stuffs

	Spectrum evalBsdf2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));
		if ((Local2::Cos(incoming) <= 0.0) || (Local2::Cos(outgoing) <= 0.0)) { return Spectrum(0.0); }

		const double exponent = mExponent->eval(vertex.mTextureCoord);
		const double halfAngle = AngleFromPhongExponent(exponent);
		const double cosHalfAngle = std::cos(halfAngle);

		const Vec2 reflect = Local2::Reflect(-incoming);
		if (Local2::Cos(reflect) <= 0.0f) { return Spectrum(0.0); }
		const double cosReflect = Math::Dot(reflect, outgoing);
		if (cosReflect < cosHalfAngle) { return Spectrum(0.0); }

		Spectrum specRefl = mReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(specRefl, Spectrum(1.0)));

		return specRefl * 0.5 / halfAngle;
	}

	double evalPdfT2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));

		// the logic here is that sample function still have a chance to generate direction under the surface.
		// we still have to evaluate the function if incoming / outgoing goes under the surface.
		if ((Local2::Cos(incoming) <= 0.0) && (Local2::Cos(outgoing) <= 0.0)) { return 0.0; }

		const double exponent = mExponent->eval(vertex.mTextureCoord);
		const double halfAngle = AngleFromPhongExponent(exponent);
		const double cosHalfAngle = std::cos(halfAngle);

		const Vec2 reflect = Local2::Reflect(-incoming);
		if (Local2::Cos(reflect) <= 0.0f) { return 0.0; }
		const double cosReflect = Math::Dot(reflect, outgoing);
		if (cosReflect < cosHalfAngle) { return 0.0; }

		return 0.5 / halfAngle;
	}

	Spectrum sampleBsdfCosByPdf2(Vec2 * outgoing, double * pdfT, const SurfaceVertex & vertex, const Vec2 & incoming, const double sample) const override
	{
		assert(Math::IsNormalized(incoming));

		if (Local2::Cos(incoming) <= 0.0) { return Spectrum(0.0); }

		const double exponent = mExponent->eval(vertex.mTextureCoord);
		const double halfAngle = AngleFromPhongExponent(exponent);
		const Vec2 reflect = Local2::Reflect(-incoming);

		const CoordFrame2 reflectBasis(reflect);
		const Vec2 out = reflectBasis.toWorld(Mapping::ArcFromSegment(sample, halfAngle));
		const double cosReflect = Math::Dot(out, reflect);
		const double cosHalfAngle = std::cos(halfAngle);

		if (outgoing) *outgoing = out;
		if (pdfT) *pdfT = 0.5 * Math::InvPi / (1.0 - cosHalfAngle);

		if (Local2::Cos(reflect) <= 0.0 || Local2::Cos(out) <= 0.0 || cosReflect <= cosHalfAngle) { return Spectrum(0.0); }
		
		Spectrum specRefl = mReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(specRefl, Spectrum(1.0)));

		return specRefl * Local2::Cos(out);
	}

	Spectrum totalEnergy2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		if (!Local2::IsFrontface(incoming) || !Local2::IsFrontface(outgoing)) return Spectrum(0.0);
		return mReflectance->eval(vertex.mTextureCoord);
	}

	// common stuffs

	BsdfType getType() const override
	{
		return BsdfType::Reflect;
	}

	shared_ptr<const Texture<Spectrum>> mReflectance;
	shared_ptr<const Texture<double>> mExponent;
};
