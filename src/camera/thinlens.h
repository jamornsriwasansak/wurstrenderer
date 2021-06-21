#pragma once

#include "common/wurst.h"

#include "common/path.h"
#include "common/util/json.h"
#include "common/camera.h"

#include "common/util/pbrt.h"

struct ThinlensCamera : public Camera
{
	static shared_ptr<ThinlensCamera> FromSimple2ParseJson(const json & json)
	{
		const Vec3 lookFrom = Util::Json::Vec3FromJson(json, "pos");
		const Vec3 lookAt = Util::Json::Vec3FromJson(json, "lookat");
		const Vec3 up = Util::Json::Vec3FromJson(json, "up");

		const double verticalFov = Util::Json::GetValue<double>(json, "fovy", 90.0) * Math::Pi / 180.0;
		if (verticalFov > Math::Pi) throw std::runtime_error("Thinlens camera doesn't support when fov > Pi radians");
		const Ivec2 resolution = Util::Json::Ivec2FromJson(json, "resolution");
		const double apertureRadius = Util::Json::GetValue<double>(json, "radius", 0.0);
		const double focalDist = Util::Json::GetValue<double>(json, "radius", 1.0);

		return make_shared<ThinlensCamera>(lookFrom, lookAt, up, verticalFov, resolution, apertureRadius, focalDist);
	}

	static shared_ptr<ThinlensCamera> FromPbrt(const pbrt::Camera & camera, const Ivec2 & resolution)
	{
		const Vec3 lookFrom = UtilPbrt::Vec3FromPbrtVec3f(camera.frame.p);
		const Vec3 lookAt = UtilPbrt::Vec3FromPbrtVec3f(camera.simplified.screen_center);
		const Vec3 up = UtilPbrt::Vec3FromPbrtVec3f(camera.simplified.screen_dv);
		const double fov = double(camera.fov) * Math::Pi / 180.0;
		const double apertureRadius = double(camera.lensRadius);
		const double focalDist = double(camera.focalDistance);

		return make_shared<ThinlensCamera>(lookFrom, lookAt, up, fov, resolution, apertureRadius, focalDist);
	}

	ThinlensCamera(const Vec3 &lookFrom, const Vec3 &lookAt, const Vec3 &up, const double verticalFov, const Ivec2 & resolution,
				   const double apertureRadius, const double focalDist, const double physicalFilmSizeByResolution = 1.0) :
		mFocalDistance(focalDist),
		mApertureRadius(apertureRadius),
		mLensArea(Math::Pi * apertureRadius * apertureRadius),
		Camera(lookFrom, verticalFov, resolution, Vec2(resolution) * physicalFilmSizeByResolution)
	{
		// basis for camera
		mCoordFrame.z = Math::Normalize(lookAt - lookFrom);
		mCoordFrame.x = Math::Normalize(Math::Cross(up, mCoordFrame.z));
		mCoordFrame.y = Math::Normalize(Math::Cross(mCoordFrame.z, mCoordFrame.x));

		// mPerFilmNormFactor = mPerPixelNormFactor * numPixels
		mPerPixelNormFactor = mDistanceLensToFilm * mDistanceLensToFilm / (mPhysicalFilmSize[0] * mPhysicalFilmSize[1]);
		mPerFilmNormFactor = mDistanceLensToFilm * mDistanceLensToFilm / (physicalFilmSizeByResolution * physicalFilmSizeByResolution);
	}

	Spectrum evalWe0(const CameraVertex &) const override
	{
		double we0 = (mApertureRadius == 0.0) ? 1.0 : 1.0 / mLensArea;
		return Spectrum(we0);
	}

	Spectrum evalWe1(Vec2 * raster, const CameraVertex & cv, const Vec3 & w) const override
	{
		const Vec3 position = mCoordFrame.toLocal(cv.mPosition - mOrigin);
		const Vec3 direction = mCoordFrame.toLocal(w);
		const double cosTheta = direction[2];

		if (Math::Dot(w, mCoordFrame.z) < 0.0) return Spectrum(0.0);

		const Vec3 focalPlanePosition = direction / cosTheta * mFocalDistance + position;
		const Vec3 filmPosition = -focalPlanePosition / mFocalDistance * mDistanceLensToFilm;
		const Vec2 rasterPosition = getFilmRasterPosition(filmPosition) * Vec2(mFilm->mResolution);
		if (!Bound2::IsInside(mFilm->mRasterBound, rasterPosition)) return Spectrum(0.0);
		if (raster) *raster = rasterPosition;

		const double cosTheta4 = Math::Square(Math::Square(cosTheta));

		return Spectrum(mPerFilmNormFactor / cosTheta4);
	}

	double evalWe0PdfA(const CameraVertex & cv) const override
	{
		return (mApertureRadius == 0.0) ? 1.0 : 1.0 / mLensArea;
	}

	double evalWe1PdfW(const CameraVertex & cv, const Vec3 & w) const override
	{
		const Vec3 position = mCoordFrame.toLocal(cv.mPosition - mOrigin);
		const Vec3 direction = mCoordFrame.toLocal(w);
		const double cosTheta = direction[2];

		#ifndef NDEBUG
            assert(Math::Dot(w, mCoordFrame.z) >= 0.0);
			const Vec3 focalPlanePosition = direction / cosTheta * mFocalDistance + position;
			const Vec3 filmPosition = -focalPlanePosition / mFocalDistance * mDistanceLensToFilm;
			const Vec2 rasterPosition = getFilmRasterPosition(filmPosition) * Vec2(mFilm->mResolution);
			assert(Bound2::IsInside(mFilm->mRasterBound, rasterPosition));
		#endif

		return mPerFilmNormFactor / (cosTheta * cosTheta * cosTheta);
	}

	double evalPerPixelWe1PdfW(const CameraVertex & cv, const Vec3 & w) const override
	{
		const Vec3 position = mCoordFrame.toLocal(cv.mPosition - mOrigin);
		const Vec3 direction = mCoordFrame.toLocal(w);
		const double cosTheta = direction[2];

		#ifndef NDEBUG
            assert(Math::Dot(w, mCoordFrame.z) >= 0.0);
			const Vec3 focalPlanePosition = direction / cosTheta * mFocalDistance + position;
			const Vec3 filmPosition = -focalPlanePosition / mFocalDistance * mDistanceLensToFilm;
			const Vec2 rasterPosition = getFilmRasterPosition(filmPosition) * Vec2(mFilm->mResolution);
			assert(Bound2::IsInside(mFilm->mRasterBound, rasterPosition));
		#endif

		return mPerPixelNormFactor / (cosTheta * cosTheta * cosTheta);
	}

	Spectrum sampleWe0ByPdf(CameraVertex * cv, double * pdfA, const Vec2 & sample) const override
	{
		// sample lens position and compute lensArea
		const Vec2 lensPoint = Mapping::DiskFromSquare(sample) * mApertureRadius;
		const Vec3 lensPosition = Vec3(lensPoint[0], lensPoint[1], 0.0);

		if (cv)
		{
			cv->mCamera = this;
			cv->mGeometryNormal = mCoordFrame.z;
			cv->mPosition = mCoordFrame.toWorld(lensPosition) + mOrigin;
		}
		if (pdfA) *pdfA = (mApertureRadius == 0.0) ? 1.0 : 1.0 / mLensArea;

		return Spectrum(1.0);
	}

	Spectrum sampleWe1CosByPdf(Vec3 * direction, Vec2 * raster, double * pdfW, const CameraVertex & cv, const Vec2 & sample) const override
	{
		// sample film position
		const Vec3 camPos = mCoordFrame.toLocal(cv.mPosition - mOrigin);
		const Vec3 filmPosition = getFilmPhysicalPosition(sample);
		if (raster) *raster = sample * Vec2(mFilm->mResolution);

		// find where ray that pass through origin intersect with focal plane
		const Vec3 focalPlanePosition = mFocalDistance * filmPosition / filmPosition[2];

		// compute outgoing direction from origin to focalPlanePosition
		const Vec3 dir = Math::Normalize(focalPlanePosition - camPos);
		const double cosTheta = dir[2];

		// transform from camera space to world space
		if (direction) *direction = mCoordFrame.toWorld(dir);
		if (pdfW) *pdfW = mPerFilmNormFactor / (cosTheta * cosTheta * cosTheta);

		return Spectrum(static_cast<double>(Math::Volume(mFilm->mResolution)));
	}

	Spectrum samplePerPixelWe1CosByPdf(Vec3 * direction, Vec2 * raster, double * pdfW, const CameraVertex & cv, const Ivec2 & pixel, const Vec2 & sample) const override
	{
		const Vec2 newSample = (Vec2(pixel) + sample) * mFilm->mInvResolution;
		sampleWe1CosByPdf(direction, raster, pdfW, cv, newSample);
		if (pdfW) *pdfW = *pdfW * Math::Volume(mFilm->mInvResolution);

		return Spectrum(1.0);
	}

	double mPerFilmNormFactor;
	double mPerPixelNormFactor;
	double mLensArea;
	double mFocalDistance;
	double mApertureRadius;
};
