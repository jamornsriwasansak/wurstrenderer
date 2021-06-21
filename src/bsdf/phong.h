#pragma once

#include "common/wurst.h"

#include "common/bsdf.h"
#include "common/path.h"
#include "common/texture.h"

// refers to Using the Modified Phong Reflectance Model for PBR by Eric P. Lafortune et al.
// http://mathinfo.univ-reims.fr/IMG/pdf/Using_the_modified_Phong_reflectance_model_for_Physically_based_rendering_-_Lafortune.pdf
// for pdf proof refers to http://www.igorsklyar.com/system/documents/papers/4/fiscourse.comp.pdf page 15
struct PhongBsdf : public Bsdf
{
	PhongBsdf(const shared_ptr<const Texture<Spectrum>> & specularReflectance, const shared_ptr<const Texture<double>> & specularLobe) :
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

		const Vec3 reflect = Local3::Reflect(-incoming);
		if (Local3::Cos(reflect) <= 0.0f) { return Spectrum(0.0); }
		const double cosReflect = Math::Dot(reflect, outgoing);
		if (cosReflect <= 0.0) { return Spectrum(0.0); }

		const double exponent = mExponent->eval(vertex.mTextureCoord);
		Spectrum specRefl = mReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(specRefl, Spectrum(1.0)));

		return specRefl * ((exponent + 2.0) * std::pow(cosReflect, exponent) * Math::InvPi * 0.5);
	}

	double evalPdfW(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assert(Math::IsNormalized(incoming));
		assert(Math::IsNormalized(outgoing));

		// the logic here is that sample function still have a chance to generate direction under the surface.
		// we still have to evaluate the function if incoming / outgoing goes under the surface.
		if ((Local3::Cos(incoming) <= 0.0) && (Local3::Cos(outgoing) <= 0.0)) { return 0.0; }

		const Vec3 reflect = Local3::Reflect(-incoming);
		const double cosReflect = Math::Dot(reflect, outgoing);
		if (cosReflect <= 0.0) { return 0.0; }

		const double exponent = mExponent->eval(vertex.mTextureCoord);
		return (exponent + 1.0) * std::pow(cosReflect, exponent) * Math::InvPi * 0.5;
	}

	BsdfType getType() const override
	{
		return BsdfType::Reflect;
	}

	Spectrum sampleBsdfCosByPdf(Vec3 * outgoing, double * pdfW, const SurfaceVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const override
	{
		assert(Math::IsNormalized(incoming));

		if (Local3::Cos(incoming) <= 0.0) { return Spectrum(0.0); }

		const double exponent = mExponent->eval(vertex.mTextureCoord);
		const Vec3 reflect = Local3::Reflect(-incoming);

		const double cosTheta = std::pow(sample[0], 1.0f / (exponent + 1.0f));
		const double sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);
		const double phi = 2.0f * Math::Pi * sample[1];
		const double cosPhi = std::cos(phi);
		const double sinPhi = std::sin(phi);
		const Vec3 localOut(sinTheta * cosPhi, cosTheta, sinTheta * sinPhi);

		const CoordFrame3 reflectBasis(reflect);
		const Vec3 out = reflectBasis.toWorld(localOut);
		const double cosReflect = Math::Dot(out, reflect);

		if (outgoing) *outgoing = out;
		if (pdfW) *pdfW = (exponent + 1.0f) * 0.5f * Math::InvPi * std::pow(cosReflect, exponent);

		if (Local3::Cos(reflect) <= 0.0 || Local3::Cos(out) <= 0.0 || cosReflect <= 0.0) { return Spectrum(0.0); }

		Spectrum specRefl = mReflectance->eval(vertex.mTextureCoord);
		assert(!Math::AnyGreater(specRefl, Spectrum(1.0)));

		// f = (n + 2) / (2 * pi) * cos^n(angle between reflect vec and sampled wi) * cosine * speccolor
		// pdf = (n + 1) / (2 * pi) * cos^n 
		return specRefl * ((exponent + 2.0f) / (exponent + 1.0f) * Local3::Cos(out));
	}

	Spectrum totalEnergy(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		if (!Local3::IsFrontface(incoming) || !Local3::IsFrontface(outgoing)) return Spectrum(0.0);
		return mReflectance->eval(vertex.mTextureCoord);
	}

	shared_ptr<const Texture<Spectrum>> mReflectance;
	shared_ptr<const Texture<double>> mExponent;
};