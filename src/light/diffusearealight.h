#pragma once

#include "common/wurst.h"

#include "common/light.h"
#include "common/geometry.h"
#include "common/path.h"
#include "common/texture.h"
#include "texture/constanttexture.h"

struct DiffuseAreaLight : public Light
{
	DiffuseAreaLight(Spectrum intensity) : mIntensity(make_shared<const ConstantTexture<Spectrum>>(intensity))
	{
	}

	DiffuseAreaLight(shared_ptr<const Texture<Spectrum>> intensity) : mIntensity(intensity)
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
		return std::max(Local3::Cos(w), 0.0) * Math::InvPi;
	}

	Spectrum evalLe0(const LightVertex & lv) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		return mIntensity->eval(lv.mTextureCoord) * Math::Pi;
	}

	Spectrum evalLe1(const LightVertex & lv, const Vec3 & w) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		if (Local3::Cos(w) < 0.0) return Spectrum(0.0);
		return Spectrum(Math::InvPi);
	}

	Spectrum sampleLe0ByPdf(LightVertex * lv, double * pdfA, const Vec2 & sample) const override
	{
		assert(lv->mGeometry->mAreaLight.get() == this);

		// sample position
		double pdf;
		SimpleVertex pointInfo = lv->mGeometry->sample(&pdf, sample);

		// populate vertex info
		lv->mCoordFrame = CoordFrame3(pointInfo.mGeometryNormal);
		lv->mTextureCoord = pointInfo.mTextureCoord;
		lv->mGeometryNormal = pointInfo.mGeometryNormal;
		lv->mLight = this;
		lv->mPosition = pointInfo.mPosition;

		// assign pdfA
		if (pdfA) *pdfA = pdf;

		return Spectrum(mIntensity->eval(lv->mTextureCoord) * Math::Pi / pdf);
	}

	Spectrum sampleLe1CosByPdf(Vec3 * outgoing, double * pdfW, const LightVertex & lv, const Vec2 & sample) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		Vec3 out = Mapping::CosineWeightedHemisphereFromSquare(sample);
		if (outgoing) *outgoing = out;
		if (pdfW) *pdfW = Local3::Cos(out) * Math::InvPi;

		return Spectrum(1.0);
	}

	// 2d stuffs

	Spectrum evalLe02(const LightVertex & lv) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		return mIntensity->eval(lv.mTextureCoord) * 2.0;
	}

	Spectrum evalLe12(const LightVertex & lv, const Vec2 & t) const override
	{
		assert(lv.mGeometry->mAreaLight.get() == this);
		if (Local2::Cos(t) < 0.0) return Spectrum(0.0);
		return Spectrum(0.5);
	}

	Spectrum sampleLe0ByPdf2(LightVertex * lv, double * pdfL, const double sample) const override
	{
		assert(lv->mGeometry->mAreaLight.get() == this);
		double pdf;
		SimpleVertex2 pointInfo = lv->mGeometry->sample2(&pdf, sample);

		lv->mPosition2 = pointInfo.mPosition2;
		lv->mGeometryNormal2 = pointInfo.mGeometryNormal2;
		lv->mTextureCoord = pointInfo.mTextureCoord;
		lv->mCoordFrame2 = CoordFrame2(pointInfo.mGeometryNormal2);
		lv->mLight = this;

		// assign pdf 
		if (pdfL) *pdfL = pdf;

		return mIntensity->eval(lv->mTextureCoord) * 2.0 / pdf;
	}

	shared_ptr<const Texture<Spectrum>> mIntensity = nullptr;
};
