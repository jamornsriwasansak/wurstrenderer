#pragma once


#include "common/wurst.h"

#include "common/camera2.h"
#include "common/path.h"

struct FlatCamera2 : public Camera2
{
	FlatCamera2(const Vec2 & start, const Vec2 & end, const int & resolution):
		mStart(start),
		mEnd(end),
		Camera2(Ivec2(resolution, 1))
	{
		Vec2 p = mEnd - mStart;
		mNormal = Math::Normalize(Vec2(p[1], p[0]));
	}

	Bound2 bound2() const override
	{
		Bound2 result;
		result = Bound2::Union(result, mStart);
		result = Bound2::Union(result, mEnd);
		return result;
	}

	Spectrum samplePerPixelWe0ByPdf(CameraVertex * cv, double * pdfL, const Ivec2 & pixel, const double sample) const override
	{
		assert(pixel[1] == 0);
		if (cv)
		{
			const double t = ((double)(pixel[0]) + sample) * mFilm->mInvResolution[0];
			cv->mCoordFrame2 = CoordFrame2(mNormal);
			cv->mCamera2 = this;
			cv->mGeometryNormal2 = mNormal;
			cv->mPosition2 = (mEnd - mStart) * t + mStart;
		}
		if (pdfL) *pdfL = mFilm->mInvResolution[0] / Math::Length(mEnd - mStart);
		return Spectrum(1.0);
	}

	Spectrum sampleWe0ByPdf(CameraVertex * cv, Vec2 * raster, double * pdfL, const double sample) const override
	{
		if (raster) *raster = Vec2(sample * (double)mFilm->mResolution[0], 0.5);
		if (cv)
		{
			cv->mCoordFrame2 = CoordFrame2(mNormal);
			cv->mCamera2 = this;
			cv->mGeometryNormal2 = mNormal;
			cv->mPosition2 = (mEnd - mStart) * sample + mStart;
		}
		if (pdfL) *pdfL = Math::Length(mEnd - mStart);
		return Spectrum((double)mFilm->mResolution[0]);
	}

	Spectrum sampleWe1CosByPdf(Vec2 * direction, double * pdfW, const CameraVertex & cv, const double sample) const override
	{
		if (direction) *direction = cv.mCoordFrame2.toWorld(Mapping::CosineWeightedHemicircleFromSegment(sample));
		if (pdfW) *pdfW = Math::InvPi;
		return Spectrum(1.0);
	}

	std::vector<std::pair<Vec2, Vec2>> segments() const override
	{
		return std::vector<std::pair<Vec2, Vec2>>(1, make_pair(mStart, mEnd));
	}

	Vec2 mNormal;
	Vec2 mStart;
	Vec2 mEnd;
};
