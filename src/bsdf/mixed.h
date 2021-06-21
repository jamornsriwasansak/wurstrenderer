#pragma once

#include <algorithm>
#include "common/bsdf.h"

struct MixedBsdf : public Bsdf
{
	MixedBsdf(std::vector<shared_ptr<const Bsdf>> & bsdfs): mBsdfs(bsdfs)
	{
	}

	// 3d stuffs

	inline void assertTotalEnergyNotExceed(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const 
	{
		Spectrum energy = totalEnergy(vertex, incoming, outgoing);
		for (Uint i = 0; i < Spectrum::NumElements; i++) { assert(energy[i] <= 1.0); }
	}

	Spectrum evalBsdf(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assertTotalEnergyNotExceed(vertex, incoming, outgoing);
		Spectrum result(0.0);
		for (const shared_ptr<const Bsdf> & bsdf : mBsdfs)
		{
			result += bsdf->evalBsdf(vertex, incoming, outgoing);
		}
		return result;
	}

	double evalPdfW(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		assertTotalEnergyNotExceed(vertex, incoming, outgoing);
		double result = 0.0;
		for (const shared_ptr<const Bsdf> & bsdf : mBsdfs)
		{
			result += bsdf->evalPdfW(vertex, incoming, outgoing);
		}
		return result / static_cast<double>(mBsdfs.size());
	}

	Spectrum sampleBsdfCosByPdf(Vec3 * outgoing, double * pdfW, const SurfaceVertex & vertex, const Vec3 & incoming, const Vec2 & sample) const override
	{
		Uint numComps = mBsdfs.size();
		Uint componentIndex = std::min((Uint)std::floor(sample[0] * numComps), numComps - 1);
		double remappedSample0 = std::min(sample[0] * numComps - componentIndex, 1.0 - Math::SmallValue);

		Vec3 out;
		Spectrum sampledBsdfCosByPdf = mBsdfs[componentIndex]->sampleBsdfCosByPdf(&out, nullptr, vertex, incoming, Vec2(remappedSample0, sample[1]));
		if (outgoing) *outgoing = out;

		assertTotalEnergyNotExceed(vertex, incoming, out);

		Spectrum bsdf = evalBsdf(vertex, incoming, out);
		double pdf = evalPdfW(vertex, incoming, out);
		if (pdfW) *pdfW = pdf;
		if (pdf == 0.0) { return Spectrum(0.0); }

		return bsdf / pdf * Local3::Cos(out);
	}

	Spectrum totalEnergy(const SurfaceVertex & vertex, const Vec3 & incoming, const Vec3 & outgoing) const override
	{
		Spectrum energy(0.0);
		for (const shared_ptr<const Bsdf> & bsdf : mBsdfs) { energy += bsdf->totalEnergy(vertex, incoming, outgoing); }
		return energy;
	}

	// 2d stuffs

	inline void assertTotalEnergyNotExceed2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const 
	{
		Spectrum energy = totalEnergy2(vertex, incoming, outgoing);
		for (Uint i = 0; i < Spectrum::NumElements; i++) { assert(energy[i] <= 1.0); }
	}

	Spectrum evalBsdf2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		assertTotalEnergyNotExceed2(vertex, incoming, outgoing);
		Spectrum result(0.0);
		for (const shared_ptr<const Bsdf> & bsdf : mBsdfs)
		{
			result += bsdf->evalBsdf2(vertex, incoming, outgoing);
		}
		return result;
	}

	double evalPdfT2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		assertTotalEnergyNotExceed2(vertex, incoming, outgoing);
		double result = 0.0;
		for (const shared_ptr<const Bsdf> & bsdf : mBsdfs)
		{
			result += bsdf->evalPdfT2(vertex, incoming, outgoing);
		}
		return result / static_cast<double>(mBsdfs.size());
	}

	Spectrum sampleBsdfCosByPdf2(Vec2 * outgoing, double * pdfT, const SurfaceVertex & vertex, const Vec2 & incoming, const double sample) const override
	{
		Uint numComps = mBsdfs.size();
		Uint componentIndex = std::min((Uint)std::floor(sample * numComps), numComps - 1);
		double remappedSample0 = std::min(sample * numComps - componentIndex, 1.0 - Math::SmallValue);

		Vec2 out;
		Spectrum sampledBsdfCosByPdf = mBsdfs[componentIndex]->sampleBsdfCosByPdf2(&out, nullptr, vertex, incoming, remappedSample0);
		if (outgoing) *outgoing = out;

		assertTotalEnergyNotExceed2(vertex, incoming, out);

		Spectrum bsdf = evalBsdf2(vertex, incoming, out);
		double pdf = evalPdfT2(vertex, incoming, out);
		if (pdfT) *pdfT = pdf;
		if (pdf == 0.0) { return Spectrum(0.0); }

		return bsdf / pdf * Local2::Cos(out);
	}

	Spectrum totalEnergy2(const SurfaceVertex & vertex, const Vec2 & incoming, const Vec2 & outgoing) const override
	{
		Spectrum energy(0.0);
		for (const shared_ptr<const Bsdf> & bsdf : mBsdfs) { energy += bsdf->totalEnergy2(vertex, incoming, outgoing); }
		return energy;
	}

	// common stuffs

	BsdfType getType() const override
	{
		BsdfType result = BsdfType(0);
		for (const shared_ptr<const Bsdf> & bsdf : mBsdfs) {
			result = result | bsdf->getType();
		}
		return result;
	}

	std::vector<shared_ptr<const Bsdf>> mBsdfs;
};
