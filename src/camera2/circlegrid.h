#pragma once

#include "common/wurst.h"

#include "common/camera2.h"
#include "common/path.h"

struct CircleGridCamera2 : public Camera2
{
	static shared_ptr<CircleGridCamera2> FromSimple2ParseJson(const json & json)
	{
		const Vec2 center = Util::Json::Vec2FromJson(json, "pos");
		const double scale = Util::Json::GetValue<double>(json, "scale", 1.0f);
		const double planeZ = Util::Json::GetValue<double>(json, "plane_z", 0.0f);
		const Ivec2 resolution = Util::Json::Ivec2FromJson(json, "resolution");

		return make_shared<CircleGridCamera2>(center, scale, planeZ, resolution);
	}

	CircleGridCamera2(const Vec2 & cameraCenter, const double scale, const double planeZ, const Ivec2 & resolution):
		Camera2(resolution),
		mCameraCenter(cameraCenter),
		mScale(scale)
	{
		mPixelSize = mFilm->mInvResolution[1] * mScale;
		const Vec2 sceneSize = mPixelSize * Vec2(resolution);
		mRenderArea = Bound2(cameraCenter - 0.5 * sceneSize, cameraCenter + 0.5 * sceneSize);
	}

	Bound2 bound2() const override
	{
		return mRenderArea;
	}

	Spectrum samplePerPixelWe0ByPdf(CameraVertex * cv, double * pdfL, const Ivec2 & pixel, const double sample) const override
	{
		const double radius = mPixelSize * 0.5;
		const Vec2 position = mPixelSize * (Vec2(pixel) + Vec2(0.5));
		const Vec2 diskPos = radius * Mapping::CircleFromSegment(sample);
		if (cv)
		{
			cv->mPosition2 = position + diskPos + mRenderArea.pMin;
			cv->mGeometryNormal2 = Vec2(0, 1);
			cv->mCamera2 = this;
		}
		if (pdfL) *pdfL = 0.5 * Math::InvPi / radius;
		return Spectrum(1.0);
	}

	Spectrum sampleWe1CosByPdf(Vec2 * direction, double * pdfW, const CameraVertex & cv, const double sample) const override
	{
		if (direction) *direction = Mapping::CircleFromSegment(sample);
		if (pdfW) *pdfW = 0.5 * Math::InvPi;
		return Spectrum(1.0);
	}

	std::vector<std::pair<Vec2, Vec2>> segments() const override
	{
		Vec2 v1(mRenderArea.pMin);
		Vec2 v2(mRenderArea.pMin[0], mRenderArea.pMax[1]);
		Vec2 v3(mRenderArea.pMin[1], mRenderArea.pMax[0]);
		Vec2 v4(mRenderArea.pMax);

		std::vector<std::pair<Vec2, Vec2>> results;
		results.push_back(make_pair(v1, v2));
		results.push_back(make_pair(v2, v4));
		results.push_back(make_pair(v1, v3));
		results.push_back(make_pair(v3, v4));

		return results;
	}

	Vec2 mCameraCenter;
	double mPixelSize;
	Bound2 mRenderArea;
	double mScale;
};
