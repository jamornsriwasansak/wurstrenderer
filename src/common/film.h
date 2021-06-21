#pragma once

#include "common/wurst.h"
#include "common/floatimage.h"
#include "common/util/enum.h"

enum class FilmRequest : uint8_t
{
	BaseBuffer = 1 << 0,
	AtomicBuffer = 1 << 1
};

DECLARE_ENUM_OPERATORS(FilmRequest);

// currently support only box filter
struct Film
{
	using AtomicSpectrum = std::array<std::atomic<double>, Spectrum::NumElements>;

	Film(const Ivec2 & size):
		mBaseScalingFactor(1.0),
		mAtomicScalingFactor(1.0),
		mResolution(size),
		mInvResolution(Vec2(1.0) / Vec2(size)),
		mPixelBound(Ivec2(0, 0), Ivec2(size)),
		mRasterBound(Vec2(0, 0), Vec2(static_cast<double>(size[0]), static_cast<double>(size[1])))
	{
	}

	Film(const Ivec2 & size, const FilmRequest & req): Film(size)
	{
		request(req);
	}

	void addSample(const Ivec2 & pixel, const Spectrum & s)
	{
		assert(mHasBaseBuffer);
		assert(Ibound2::IsInside(mPixelBound, pixel));
		getBasePixel(pixel) += s;
	}

	void addSample(const Vec2 & raster, const Spectrum & s)
	{
		assert(mHasBaseBuffer);
		assert(Bound2::IsInside(mRasterBound, raster));
		Ivec2 p(raster);
		addSample(p, s);
	}

	void atomicAddSample(const Ivec2 & pixel, const Spectrum & s)
	{
		assert(mHasAtomicBuffer);
		assert(Ibound2::IsInside(mPixelBound, pixel));
		AtomicSpectrum & as = getAtomicPixel(pixel);
		for (Uint iChannel = 0; iChannel < Spectrum::NumElements; iChannel++)
		{
			double prev = as[iChannel].load();
			while (!as[iChannel].compare_exchange_weak(prev, prev + s[iChannel]));
		}
	}

	void atomicAddSample(const Vec2 & raster, const Spectrum & s)
	{
		assert(mHasAtomicBuffer);
		assert(Bound2::IsInside(mRasterBound, raster));
		Ivec2 p(raster);
		atomicAddSample(p, s);
	}

	void clear()
	{
		for (Spectrum & s : mBaseBuffer) { s = Spectrum(0.0); }
		for (AtomicSpectrum & as : mAtomicBuffer)
			for (Uint iChannel = 0; iChannel < Spectrum::NumElements; iChannel++)
				as[iChannel].store(0.0);
	}

	Fimage<Spectrum> getImage() const
	{
		Fimage<Spectrum> result(mResolution);
		const double atomicScalingFactor = mAtomicScalingFactor.load();
		result.forEachPixel([&](const Ivec2 & pixel, Spectrum * color)
		{
			Spectrum c(0.0);
			if (mHasBaseBuffer) c += getBasePixel(pixel) * mBaseScalingFactor;
			if (mHasAtomicBuffer) c += readAtomicPixel(pixel) * atomicScalingFactor;
			*color = c;
		});
		return result;
	}

	AtomicSpectrum & getAtomicPixel(const Ivec2 & pixel)
	{
		assert(mHasAtomicBuffer);
		assert(Ibound2::IsInside(mPixelBound, pixel));
		return mAtomicBuffer[pixel[1] * mResolution[0] + pixel[0]];
	}

	Spectrum readAtomicPixel(const Ivec2 & pixel) const
	{
		assert(mHasAtomicBuffer);
		assert(Ibound2::IsInside(mPixelBound, pixel));
		const AtomicSpectrum & as = mAtomicBuffer[pixel[1] * mResolution[0] + pixel[0]];
		Spectrum result;
		for (Uint i = 0; i < Spectrum::NumElements; i++) { result[i] = as[i].load(); }
		return result;
	}

	Spectrum & getBasePixel(const Ivec2 & pixel)
	{
		assert(mHasBaseBuffer);
		assert(Ibound2::IsInside(mPixelBound, pixel));
		return mBaseBuffer[pixel[1] * mResolution[0] + pixel[0]];
	}

	Spectrum getBasePixel(const Ivec2 & pixel) const
	{
		assert(mHasBaseBuffer);
		assert(Ibound2::IsInside(mPixelBound, pixel));
		return mBaseBuffer[pixel[1] * mResolution[0] + pixel[0]];
	}

	Spectrum getPixel(const Ivec2 & pixel) const
	{
		assert(Ibound2::IsInside(mPixelBound, pixel));
		Spectrum result(0.0);
		if (mHasBaseBuffer) result += getBasePixel(pixel) * mBaseScalingFactor;
		if (mHasAtomicBuffer) result += readAtomicPixel(pixel) * mAtomicScalingFactor.load();
		return result;
	}

	void request(const FilmRequest request)
	{
		if ((request & FilmRequest::BaseBuffer) != 0) requestBaseBuffer();
		if ((request & FilmRequest::AtomicBuffer) != 0) requestAtomicBuffer();
	}

	void requestBaseBuffer()
	{
		mBaseBuffer = std::vector<Spectrum>(mResolution[0] * mResolution[1]);
		mHasBaseBuffer = true;
		resetBaseBuffer();
	}

	void requestAtomicBuffer()
	{
		mAtomicBuffer = std::vector<AtomicSpectrum>(mResolution[0] * mResolution[1]);
		mHasAtomicBuffer = true;
		resetAtomicBuffer();
	}

	void setBasePixel(const Ivec2 & pixel, const Spectrum & s)
	{
		assert(mHasBaseBuffer);
		assert(Ibound2::IsInside(mPixelBound, pixel));
		mBaseBuffer[pixel[1] * mResolution[0] + pixel[0]] = s;
	}

	void setAtomicPixel(const Ivec2 & pixel, const Spectrum & s)
	{
		assert(mHasAtomicBuffer);
		assert(Ibound2::IsInside(mPixelBound, pixel));
		AtomicSpectrum & as = mAtomicBuffer[pixel[1] * mResolution[0] + pixel[0]];
		for (Uint i = 0; i < Spectrum::NumElements; i++) { as[i].store(s[i]); }
	}

	void resetAtomicBuffer()
	{
		for (AtomicSpectrum & as : mAtomicBuffer) for (std::atomic<double> & s : as) { s.store(0.0); }
	}

	void resetBaseBuffer()
	{
		for (Spectrum & s : mBaseBuffer) { s = Spectrum(0.0); }
	}

	std::string mName = "";
	bool mHasBaseBuffer = false;
	bool mHasAtomicBuffer = false;
	std::vector<Spectrum> mBaseBuffer;
	std::vector<AtomicSpectrum> mAtomicBuffer;
	Ivec2 mResolution;
	Bound2 mRasterBound;
	Ibound2 mPixelBound;
	Vec2 mInvResolution;
	double mBaseScalingFactor;
	std::atomic<double> mAtomicScalingFactor;
};
