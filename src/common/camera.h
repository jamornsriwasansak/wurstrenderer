#pragma once

#include "common/wurst.h"

#include "common/film.h"

struct Camera
{
	Camera(const Vec3 & origin, const double verticalFov, const Ivec2 & resolution, const Vec2 & physicalFilmSize) :
		mOrigin(origin),
		mVerticalFov(verticalFov),
		mPhysicalFilmSize(physicalFilmSize),
		mFilm(make_shared<Film>(resolution))
	{
		mFilm->clear();
		mDistanceLensToFilm = mPhysicalFilmSize[1] * 0.5 / std::tan(verticalFov * 0.5);
	}

	// film physical position is already filpped
	Vec3 getFilmPhysicalPosition(const Vec2 & uv) const
	{
		const Vec2 ndc = uv - Vec2(0.5);
		const Vec2 xy = ndc * mPhysicalFilmSize;
		return Vec3(xy[0], xy[1], -mDistanceLensToFilm);
	}

	Vec2 getFilmRasterPosition(const Vec3 & position) const
	{
		assert(Math::IsApprox(position[2], -mDistanceLensToFilm));
		const Vec2 xy(position[0], position[1]);
		const Vec2 ndc = xy / mPhysicalFilmSize;
		return ndc + Vec2(0.5);
	}

	virtual Spectrum evalWe0(const CameraVertex & cv) const = 0;

	virtual Spectrum evalWe1(Vec2 * raster, const CameraVertex & cv, const Vec3 & w) const = 0;

	virtual double evalWe0PdfA(const CameraVertex & cv) const = 0;

	virtual double evalWe1PdfW(const CameraVertex & cv, const Vec3 & w) const = 0;

	virtual double evalPerPixelWe1PdfW(const CameraVertex & cv, const Vec3 & w) const = 0;

	virtual Spectrum sampleWe0ByPdf(CameraVertex * cv, double * pdfA, const Vec2 & sample) const = 0;

	virtual Spectrum sampleWe1CosByPdf(Vec3 * direction, Vec2 * raster, double * pdfW, const CameraVertex & cv, const Vec2 & sample) const = 0;

	virtual Spectrum samplePerPixelWe1CosByPdf(Vec3 * direction, Vec2 * raster, double * pdfW, const CameraVertex & cv, const Ivec2 & pixel, const Vec2 & sample) const = 0;

	Vec2 mPhysicalFilmSize;
	Vec3 mOrigin;
	CoordFrame3 mCoordFrame;
	shared_ptr<Film> mFilm = nullptr;
	double mDistanceLensToFilm = 0.0;
	double mVerticalFov = 0.0;
	mutable shared_ptr<const Medium> mMedium = nullptr;
};
