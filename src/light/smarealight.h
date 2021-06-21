#pragma once

#include "common/wurst.h"

#include "common/light.h"
#include "common/geometry.h"
#include "common/path.h"
#include "common/texture.h"
#include "texture/constanttexture.h"

// same as diffuse area light 
struct SimpleMetalAreaLight : public Light
{
	SimpleMetalAreaLight(Spectrum intensity, double halfAngle):
		mIntensity(make_shared<const ConstantTexture<Spectrum>>(intensity)),
		mHalfAngle(halfAngle),
		mCosHalfAngle(std::cos(halfAngle))
	{
	}

	SimpleMetalAreaLight(shared_ptr<const Texture<Spectrum>> intensity, const double halfAngle):
		mIntensity(intensity),
		mHalfAngle(halfAngle),
		mCosHalfAngle(std::cos(halfAngle))
	{
	}
	
	double cdfWeight() const override
	{
		return Math::Length(mIntensity->average());
	}

	// 3d stuffs

	double evalPdfLe0(const LightVertex & lv) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		return lv.mGeometry->evalPdfA();
	}

	double evalPdfLe1(const LightVertex & lv, const Vec3 & w) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		if (Local3::Cos(w) < mCosHalfAngle) return 0.0;
		return 0.5 * Math::InvPi / (1.0 - mCosHalfAngle);
	}

	Spectrum evalLe0(const LightVertex & lv) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		return mIntensity->eval(lv.mTextureCoord) * Math::Pi;
	}

	Spectrum evalLe1(const LightVertex & lv, const Vec3 & w) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		if (Local3::Cos(w) < mCosHalfAngle) return Spectrum(0.0);
		return Spectrum(0.5 * Math::InvPi / (1.0 - mCosHalfAngle));
	}

	Spectrum sampleLe0ByPdf(LightVertex * lv, double * pdfA, const Vec2 & sample) const override
	{
		assert(lv->mGeometry->mAreaLight.get() == this);

		// sample position
		double pdf;
		SimpleVertex pointInfo = lv->mGeometry->sample(&pdf, sample);

		// populate vertex info
		if (lv)
		{
			lv->mCoordFrame = CoordFrame3(pointInfo.mGeometryNormal);
			lv->mGeometryNormal = pointInfo.mGeometryNormal;
			lv->mLight = this;
			lv->mPosition = pointInfo.mPosition;
		}

		// assign pdfA
		if (pdfA) *pdfA = pdf;

		return Spectrum(mIntensity->eval(lv->mTextureCoord) * Math::Pi / pdf);
	}

	Spectrum sampleLe1CosByPdf(Vec3 * outgoing, double * pdfW, const LightVertex & lv, const Vec2 & sample) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		const Vec3 out = Mapping::SolidAngleFromSquare(sample, mHalfAngle);
		if (outgoing) *outgoing = out;
		if (pdfW) *pdfW = 0.5 * Math::InvPi / (1.0 - mCosHalfAngle);
		return Spectrum(Local3::Cos(out));
	}

	// 2d stuffs

	Spectrum evalLe02(const LightVertex & lv) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		return mIntensity->eval(lv.mTextureCoord) * 2.0 * mHalfAngle;
	}
	
	Spectrum evalLe12(const LightVertex & lv, const Vec2 & w) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		if (Local2::Cos(w) < mCosHalfAngle) return Spectrum(0.0);
		return Spectrum(0.5 / mHalfAngle);
	}

	shared_ptr<const Texture<Spectrum>> mIntensity = nullptr;
	double mHalfAngle;
	double mCosHalfAngle;
};
