#pragma once

#include "vecmath.h"

#if __AVX2__ 
using RgbSpectrum = Math::VecSimd4<3, 1001>;
#else
using RgbSpectrum = Math::VecPrim<double, 3, 1001>;
#endif

struct SpectrumConvertUtil
{
	static inline double FromSrgb(const double v)
	{
		if (v <= 0.04045) return v * (1.0 / 12.92);
		return std::pow((v + 0.055) * (1.0 / 1.055), 2.4);
	}

	static RgbSpectrum SpectrumFromSrgb(const double r, const double g, const double b)
	{
		return RgbSpectrum(FromSrgb(r), FromSrgb(g), FromSrgb(b));
	}

	static RgbSpectrum SpectrumFromRgb(const double r, const double g, const double b)
	{
		return RgbSpectrum(r, g, b);
	}

	static RgbSpectrum SpectrumFromRg(const Vec2 & v)
	{
		return RgbSpectrum(v[0], v[1], 0.0);
	}

	static RgbSpectrum SpectrumFromRgb(const Vec3 & v)
	{
		return RgbSpectrum(v[0], v[1], v[2]);
	}

	static Vec3 Vec3FromSpectrum(const RgbSpectrum & spectrum)
	{
		return Vec3(spectrum[0], spectrum[1], spectrum[2]);
	}

	static std::vector<RgbSpectrum> SpectrumsFromRgbs(const std::vector<Vec3> & rgbs)
	{
		std::vector<RgbSpectrum> results(rgbs.size());
		for (Uint i = 0;i < rgbs.size();i++) { results[i] = SpectrumFromRgb(rgbs[i]); }
		return results;
	}

	static std::vector<Vec3> Vec3sFromSpectrums(const std::vector<RgbSpectrum> & spectrums)
	{
		std::vector<Vec3> results(spectrums.size());
		for (Uint i = 0;i < spectrums.size();i++) { results[i] = Vec3FromSpectrum(spectrums[i]); }
		return results;
	}
};
