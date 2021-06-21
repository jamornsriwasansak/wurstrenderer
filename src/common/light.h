#pragma once

#include "common/wurst.h"

#include "common/cdftable.h"

struct Light
{
	Light() {};

	virtual double cdfWeight() const = 0;

	// 3d

	virtual double evalPdfLe0(const LightVertex & vertex) const = 0;

	virtual double evalPdfLe1(const LightVertex & vertex, const Vec3 & outgoing) const { throw std::runtime_error("unimpl"); return 0.0; } //TODO:: implement this for all classses

	virtual Spectrum evalLe0(const LightVertex & vertex) const = 0;

	virtual Spectrum evalLe1(const LightVertex & vertex, const Vec3 & w) const = 0;

	virtual Spectrum sampleLe0ByPdf(LightVertex * lightVertex, double * pdfA, const Vec2 & sample) const = 0;

	virtual Spectrum sampleLe1CosByPdf(Vec3 * out, double * pdfW, const LightVertex & vertex, const Vec2 & sample) const = 0;

	// 2d

	virtual Spectrum evalLe02(const LightVertex & vertex) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	virtual Spectrum evalLe12(const LightVertex & vertex, const Vec2 & t) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	virtual Spectrum sampleLe0ByPdf2(LightVertex * lv, double * pdfL, const double sample) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	virtual Spectrum sampleLe1CosByPdf2(Vec2 * outgoing, double * pdfT, const LightVertex & vertex, const double sample) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	shared_ptr<const Medium> mMedium = nullptr;
};
