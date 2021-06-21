#pragma once

#include "common/wurst.h"

#include "common/cdftable.h"
#include "common/floatimage.h"
#include "common/mapping.h"
#include "common/parallel.h"
#include "common/sphericalrepresentation.h"
#include "light/envmap.h"

template <typename Type>
struct EquiRectRep : public SphRep<Type>
{
	EquiRectRep()
	{}

	EquiRectRep(const Type & s)
	{
		mImage = Fimage<Type>(1, 1);
		mImage.mWrapS = Fimage<Type>::WrapMode::Repeat;
		mImage.mWrapT = Fimage<Type>::WrapMode::Clamp;
		mImage.at(0, 0) = s;
	}

	EquiRectRep(const Fimage<Type> & fimage): mImage(fimage)
	{
		mImage.mWrapS = Fimage<Type>::WrapMode::Repeat;
		mImage.mWrapT = Fimage<Type>::WrapMode::Clamp;
		if (fimage.mSize[0] != fimage.mSize[1] * 2) throw std::runtime_error("envmap has wrong size");
	}

	EquiRectRep(const Fimage<Type> & fimage, const bool doBuildCdfTable) :
		mImage(fimage)
	{
		mImage.mWrapS = Fimage<Type>::WrapMode::Repeat;
		mImage.mWrapT = Fimage<Type>::WrapMode::Clamp;
		if (fimage.mSize[0] != fimage.mSize[1] * 2) throw std::runtime_error("envmap has wrong size");

		if (doBuildCdfTable)
		{
			std::vector<double> weights(fimage.getNumPixels());
			// divide by envmap radius seems to be good enough
			for (int y = 0; y < fimage.mSize[1]; y++)
			{
				const double v = (static_cast<double>(y) + 0.5) / static_cast<double>(fimage.mSize[1]);
				const double sine = std::sin(v * Math::Pi);
				for (int x = 0; x < fimage.mSize[0]; x++)
				{
					Type approxSpec = fimage.evalTexel(x, y);

					// pixel area on the sphere = sine * #rows / Pi * #cols / (2 Pi)
					// and we can ignore the constant part since they will eventually be normalized in cdftable
					const double unnormalizedArea = sine;
					const double approxLuminance = Math::Length(approxSpec);

					// must add small value incase of fimage with special filter to prevent the case where some pixels aren't sampled at all
					const double weight = unnormalizedArea * approxLuminance;
					weights[y * fimage.mSize[0] + x] = weight;
				}
			};
			mCdfTable = CdfTable(weights);
		}
	}

	double pdfA(const Vec3 & dir) const override
	{
		Ivec2 pixel = Ivec2(Mapping::PanoramaFromWorld(dir) * Vec2(mImage.mSize));
		int pixelIndex = mImage.mSize[0] * pixel[1] + pixel[0];
		const double sine = Local3::Sin(dir);
		if (Math::IsApprox(sine, 0.0)) return 0.0;
		return mCdfTable.mPdfs[pixelIndex] * mImage.getNumPixels() / sine * Math::InvPi * Math::InvPi * 0.5;
	}

	Type eval(const Vec3 & dir) const override
	{
		Vec2 textureUv = Mapping::PanoramaFromWorld(dir);
		return mImage.evalBilinear(textureUv);
	}

	EquiRectRep<Type> rotate(const CoordFrame3 & frame) const
	{
		Fimage<Type> rotatedImage(mImage.mSize);
		rotatedImage.forEachPixel([&](const Vec2 & uv, double * color) {
			const Vec3 direction = Mapping::WorldFromPanorama(uv);
			const Vec3 trueDirection = frame.toLocal(direction);
			*color = this->eval(trueDirection);
		});
		return EquiRectRep<Type>(rotatedImage);
	}

	Type sample(Vec3 * dir, double * pdfW, const Vec2 & sample) const override
	{
		if (mCdfTable.isEmpty())
		{
			// cdftable is not constructed use the uniform sphere sampling method.
			const Vec3 normal = Mapping::SphereFromSquare(sample);
			if (dir) *dir = normal;
			if (pdfW) *pdfW = 0.25 * Math::InvPi;
			return eval(-normal) * 4.0 * Math::Pi;
		}
		else
		{
			Vec2 sample2 = sample;
			double pdf;
			Uint pixelIndex = mCdfTable.sample(&pdf, &sample2[0], sample2[0]);
			Uint x = pixelIndex % mImage.mSize[0];
			Uint y = pixelIndex / mImage.mSize[0];

			const Vec3 direction = Mapping::WorldFromPanorama(Vec2((double)x + sample2[0], (double)y + sample2[1]) / Vec2(mImage.mSize));
			const double sine = Local3::Sin(direction);

			// this pdf consists of two parts
			// 1. randomly sample the pixel
			// 2. randomly sample inside that pixel (pixel area on the sphere = sine * (Pi / Y) * (2 * Pi / X)
			pdf = Math::IsApprox(sine, 0.0) ? 0.0 : pdf * mImage.getNumPixels() / sine * Math::InvPi * Math::InvPi * 0.5;

			if (dir) *dir = direction;
			if (pdfW) *pdfW = pdf;

			return eval(direction) / pdf;
		}
	}

	Fimage<Type> mImage;
	CdfTable mCdfTable;
};
