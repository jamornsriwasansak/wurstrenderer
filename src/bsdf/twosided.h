#pragma once

#include "common/wurst.h"

#include "common/bsdf.h"
#include "common/path.h"

struct TwoSidedBsdf : public Bsdf
{
	TwoSidedBsdf(const shared_ptr<const Bsdf> & bsdfBackface, const shared_ptr<const Bsdf> & bsdfFrontface):
		mBsdfBackface(bsdfBackface), mBsdfFrontface(bsdfFrontface)
	{
		assert(bsdfBackface->getType() == BsdfType::Reflect);
		assert(bsdfFrontface->getType() == BsdfType::Reflect);
	}

	TwoSidedBsdf(const shared_ptr<const Bsdf> & bsdf):
		mBsdfBackface(bsdf), mBsdfFrontface(bsdf)
	{
		assert(bsdf->getType() == BsdfType::Reflect);
	}

	// 3d stuffs

	Spectrum evalBsdf(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		if (Local3::IsFrontface(incoming) && Local3::IsFrontface(outgoing))
		{
			return mBsdfFrontface->evalBsdf(vertex, incoming, outgoing);
		}
		else if (!Local3::IsFrontface(incoming) && !Local3::IsFrontface(outgoing))
		{
			return mBsdfBackface->evalBsdf(vertex, Math::FlipY(incoming), Math::FlipY(outgoing));
		}
		return Spectrum(0.0);
	}

	double evalPdfW(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		if (Local3::IsFrontface(incoming) && Local3::IsFrontface(outgoing))
		{
			return mBsdfFrontface->evalPdfW(vertex, incoming, outgoing);
		}
		else if (!Local3::IsFrontface(incoming) && !Local3::IsFrontface(outgoing))
		{
			return mBsdfBackface->evalPdfW(vertex, Math::FlipY(incoming), Math::FlipY(outgoing));
		}
		return 0.0;
	}

	Spectrum sampleBsdfCosByPdf(Vec3 * outgoing, double * pdfW, const SurfaceVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const override
	{
		if (Local3::IsFrontface(incoming))
		{
			return mBsdfFrontface->sampleBsdfCosByPdf(outgoing, pdfW, vertex, incoming, sample);
		}
		else
		{
			Spectrum result = mBsdfBackface->sampleBsdfCosByPdf(outgoing, pdfW, vertex, Math::FlipY(incoming), sample);
			if (outgoing) *outgoing = Math::FlipY(*outgoing);
			return result;
		}
	}

	Spectrum totalEnergy(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		if (Local3::IsFrontface(incoming) && Local3::IsFrontface(outgoing))
		{
			return mBsdfFrontface->totalEnergy(vertex, incoming, outgoing);
		}
		else if (!Local3::IsFrontface(incoming) && !Local3::IsFrontface(outgoing))
		{
			return mBsdfBackface->totalEnergy(vertex, incoming, outgoing);
		}
		return Spectrum(0.0);
	}

	// 2d stuffs

	Spectrum evalBsdf2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		if (Local2::IsFrontface(incoming) && Local2::IsFrontface(outgoing))
		{
			return mBsdfFrontface->evalBsdf2(vertex, incoming, outgoing);
		}
		else if (!Local2::IsFrontface(incoming) && !Local2::IsFrontface(outgoing))
		{
			return mBsdfBackface->evalBsdf2(vertex, -incoming, -outgoing);
		}
		return Spectrum(0.0);
	}

	double evalPdfT2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		if (Local2::IsFrontface(incoming) && Local2::IsFrontface(outgoing))
		{
			return mBsdfFrontface->evalPdfT2(vertex, incoming, outgoing);
		}
		else if (!Local2::IsFrontface(incoming) && !Local2::IsFrontface(outgoing))
		{
			return mBsdfBackface->evalPdfT2(vertex, -incoming, -outgoing);
		}
		return 0.0;
	}

	Spectrum sampleBsdfCosByPdf2(Vec2 * outgoing, double * pdfT, const SurfaceVertex & vertex, const Vec2 & incoming, const double sample) const override
	{
		if (Local2::IsFrontface(incoming))
		{
			return mBsdfFrontface->sampleBsdfCosByPdf2(outgoing, pdfT, vertex, incoming, sample);
		}
		else
		{
			Spectrum result = mBsdfBackface->sampleBsdfCosByPdf2(outgoing, pdfT, vertex, -incoming, sample);
			if (outgoing) *outgoing = -*outgoing;
			return result;
		}
	}

	Spectrum totalEnergy2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		if (Local2::IsFrontface(incoming) && Local2::IsFrontface(outgoing))
		{
			return mBsdfFrontface->totalEnergy2(vertex, incoming, outgoing);
		}
		else if (!Local2::IsFrontface(incoming) && !Local2::IsFrontface(outgoing))
		{
			return mBsdfBackface->totalEnergy2(vertex, incoming, outgoing);
		}
		return Spectrum(0.0);
	}

	// common stuffs

	BsdfType getType() const override
	{
		return BsdfType::Reflect | BsdfType::TwoSided;
	}

	shared_ptr<const Bsdf> mBsdfBackface;
	shared_ptr<const Bsdf> mBsdfFrontface;
};
