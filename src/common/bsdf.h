#pragma once

#include "common/wurst.h"
#include "common/util/file.h"
#include "common/util/enum.h"

enum class BsdfType : uint8_t
{
	Refract = 2 << 0,
	Reflect = 2 << 1,
	TwoSided = 2 << 2
};

DECLARE_ENUM_OPERATORS(BsdfType);

struct Bsdf
{
	// 3d stuffs

	virtual Spectrum evalBsdf(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const = 0;

	virtual double evalPdfW(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const = 0;

	virtual Spectrum sampleBsdfCosByPdf(Vec3 * outgoing, double * pdfW, const SurfaceVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const = 0;

	virtual Spectrum totalEnergy(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const = 0;

	// 2d stuffs

	virtual Spectrum evalBsdf2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	virtual double evalPdfT2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const { throw std::runtime_error("unimpl"); return 0.0; }

	virtual Spectrum sampleBsdfCosByPdf2(Vec2 * outgoing, double * pdfT, const SurfaceVertex & vertex, const Vec2 & incoming, const double sample) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	virtual Spectrum totalEnergy2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	// common

	virtual BsdfType getType() const = 0;

	virtual bool hasType(const BsdfType & type) const { return (this->getType() & type) != 0; }
};