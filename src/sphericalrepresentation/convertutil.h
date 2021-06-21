#pragma once

#include "sphericalrepresentation/equirectangular.h"
#include "sphericalrepresentation/sphericalharmonics.h"
#include "sphericalrepresentation/sphericalfibonacci.h"

struct SphRepConvertUtil
{
	template <typename Type>
	static EquiRectRep<Type> EquiRectRepFrom(const SphHarRep<Type> & shrep, const Ivec2 & resolution) 
	{
		Fimage<Type> image(resolution);
		for (int iCoeff = 0; iCoeff < shrep.mCoeffs.size(); iCoeff++)
		{
			image.forEachPixel([&](const Vec2 & uv, Type * color) {
				const double shVal = Math::SphericalHarmonics::Eval(iCoeff, Mapping::WorldFromPanorama(uv));
				*color += shVal * shrep.mCoeffs[iCoeff];
			});
		}
		return EquiRectRep(image);
	}

	template <typename Type>
	static EquiRectRep<Type> EquiRectRepFrom(const SphFiboRep<Type> & sfrep, const Ivec2 & resolution)
	{
		Fimage<Type> image(resolution);
		image.forEachPixel([&](const Vec2 & uv, Type * color) {
			*color = sfrep.eval(Mapping::WorldFromPanorama(uv));
		});
		return EquiRectRep(image);
	}

	template <typename Type>
	static SphFiboRep<Type> SphFiboRepFrom(const EquiRectRep<Type> & errep, const int resolution)
	{
		SphFiboRep<Type> result(resolution);
		for (int i = 0; i < resolution; i++)
		{
			const Vec3 dir = Mapping::WorldFromSphFiboIndex(i, resolution);
			const Vec2 uv = Mapping::PanoramaFromWorld(dir);
			result.mData[i] = errep.mImage.evalBilinear(uv);
		}
		return result;
	}

	template <typename Type>
	static SphFiboRep<Type> SphFiboRepFrom(const EquiRectRep<Type> & errep, const int resolution, Sampler * sampler, const int numSamplesPerCell)
	{
		SphFiboRep<Type> result(resolution);

		// approximate that each spherical fibonacci block takes a shape of cone
		const double solidAngle = Math::Pi * 4.0 / static_cast<double>(resolution);
		const double halfConeAngle = std::acos(1.0 - solidAngle * Math::InvPi * 0.5);

		for (int i = 0; i < resolution; i++)
		{
			const Vec3 dir = Mapping::WorldFromSphFiboIndex(i, resolution);
			const CoordFrame3 basis(dir);
			Type sum(0.0);
			for (int j = 0; j < numSamplesPerCell; j++)
			{
				Vec3 sampledDir = basis.toWorld(Mapping::SolidAngleFromSquare(sampler->get2d(), halfConeAngle));
				const Vec2 uv = Mapping::PanoramaFromWorld(dir);
				sum += errep.mImage.evalBilinear(uv);
			}
			result.mData[i] = sum / static_cast<double>(numSamplesPerCell);
		}
		return result;
	}

	template <typename Type>
	static SphFiboRep<Type> SphFiboRepFrom(const SphHarRep<Type> & sphHarRep, const int resolution)
	{
		SphFiboRep<Type> result(resolution);
		for (int i = 0; i < resolution; i++)
		{
			const Vec3 dir = Mapping::WorldFromSphFiboIndex(i, resolution);
			result.mData[i] = sphHarRep.eval(dir);
		}
		return result;
	}

	template <typename Type>
	static SphHarRep<Type> SphHarRepFrom(const EquiRectRep<Type> & equirect, const int numCoeffs)
	{
		const double invNumPixels = 1.0 / static_cast<double>(equirect.mImage.getNumPixels());
		const Fimage<Type> & image = equirect.mImage;

		SphHarRep<Type> result(numCoeffs);
		double normalizationFactor = 2.0 * Math::Pi * Math::Pi * invNumPixels;
		for (int iCoeff = 0; iCoeff < (int)numCoeffs; iCoeff++)
		{
			Type coeff(0.0);
			image.forEachPixel([&](const Vec2 & uv, const Type & color) {
				const Vec2 thetaPhi = uv * Vec2(2.0 * Math::Pi, Math::Pi);
				double shVal = Math::SphericalHarmonics::Eval(iCoeff, thetaPhi[0], thetaPhi[1]);
				coeff += shVal * color * std::sin(thetaPhi[1]) * normalizationFactor;
			});
			result.mCoeffs[iCoeff] = coeff;
		}
		return result;
	}

	template <typename Type>
	static SphHarRep<Type> SphHarRepFrom(const SphFiboRep<Type> & sfrep, const int numCoeffs)
	{
		const double invNumFiboSamples = 1.0 / static_cast<double>(sfrep.mData.size());
		SphHarRep<Type> result(numCoeffs);
		for (int iCoeff = 0; iCoeff < numCoeffs; iCoeff++)
		{
			Type coeff(0.0);
			for (int iFibo = 0; iFibo < sfrep.mData.size(); iFibo++)
			{
				const Vec3 direction = Mapping::WorldFromSphFiboIndex(iFibo, sfrep.mData.size());
				const Vec2 thetaPhi = Mapping::SphericalFromWorld(direction);
				const double shVal = Math::SphericalHarmonics::Eval(iCoeff, thetaPhi[0], thetaPhi[1]);
				coeff += shVal * sfrep.mData[iFibo];
			}
			coeff = coeff * invNumFiboSamples * Math::Pi * 4.0;
			result.mCoeffs[iCoeff] = coeff;
		}
		return result;
	}
};