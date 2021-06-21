#pragma once

#include "common/wurst.h"

#include "common/film.h"
#include "common/floatimage.h"
#include "common/path.h"
#include "common/geometry.h"
#include "common/scene.h"

#include "camera2/flat.h"

struct Visualizer2
{
	Visualizer2(const shared_ptr<const Camera2> & camera2, const shared_ptr<const Scene> & scene, const Ivec2 & resolution):
		mCamera2(camera2),
		mNumDrawnPaths(0),
		mScene(scene),
		mFilm(make_shared<Film>(resolution))
	{
		mBound2 = Bound2::Union(camera2->bound2(), scene->bound2());
		mBound2 = Bound2::Resize(mBound2, 1.1);
		mFilm->requestBaseBuffer();
		mFilm->requestAtomicBuffer();
		drawCamera2();
		drawScene();
	}

	void drawCamera2(const Spectrum s = SpectrumConvertUtil::SpectrumFromRgb(1.0, 0.0, 0.0))
	{
		std::vector<std::pair<Vec2, Vec2>> segments = mCamera2->segments();
		for (const std::pair<Vec2, Vec2> & segment : segments)
		{
			this->drawLine(Ivec2(rasterFromWorld2(segment.first)), Ivec2(rasterFromWorld2(segment.second)), false, s);
		}
	}

	void drawLine(const Ivec2 & a, const Ivec2 & b, const bool atomic, const Spectrum & s = Spectrum(1.0))
	{
		int x0 = (int)a[0], x1 = (int)b[0], y0 = (int)a[1], y1 = (int)b[1];
		int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
		int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
		int err = (dx > dy ? dx : -dy) / 2, e2;
		while (true)
		{
			if (Ibound2::IsInside(Ibound2(Ivec2(0, 0), Ivec2(mFilm->mResolution)), Ivec2(x0, y0)))
			{
				if (atomic)
				{
					mFilm->atomicAddSample(Ivec2(x0, y0), s);
				}
				else
				{
					mFilm->addSample(Ivec2(x0, y0), s);
				}
			}
			if (x0 == x1 && y0 == y1) break;
			e2 = err;
			if (e2 > -dx) { err -= dy; x0 += sx; }
			if (e2 < dy) { err += dx; y0 += sy; }
		}
	}

	void drawPath(const std::vector<shared_ptr<Vertex>> & path, bool doClearPreviousPath = false)
	{
		if (doClearPreviousPath)
		{
			mFilm->resetAtomicBuffer();
			mNumDrawnPaths = 0;
		}

		const Spectrum normalColor = SpectrumConvertUtil::SpectrumFromRgb(0.0, 5.0, 0.0);
		for (const shared_ptr<Vertex> & vertex : path)
		{
			Vec2 p = vertex->mPosition2;
			Vec2 n = vertex->mGeometryNormal2;
			this->drawLine(Ivec2(rasterFromWorld2(p + n * 0.01)), Ivec2(rasterFromWorld2(p + n * 0.02)), true, normalColor);
		}

		for (Uint iVertex = 1; iVertex < path.size(); iVertex++)
		{
			Vec2 p = path[iVertex - 1]->mPosition2;
			Vec2 n = path[iVertex - 1]->mGeometryNormal2;
			Vec2 pNext = path[iVertex]->mPosition2;
			Vec2 nNext = path[iVertex]->mGeometryNormal2;
			this->drawLine(Ivec2(rasterFromWorld2(p + n * 0.001)), Ivec2(rasterFromWorld2(pNext + nNext * 0.001)), true);
		}
		mNumDrawnPaths += 1;
		mFilm->mAtomicScalingFactor.store(std::min(1.0 / static_cast<double>(mNumDrawnPaths) * 20.0, 1.0));
	}

	void drawScene()
	{
		for (const shared_ptr<const Geometry> & geometry : mScene->mGeometries)
		{
			Spectrum color = (geometry->mAreaLight) ? SpectrumConvertUtil::SpectrumFromRgb(1.0, 1.0, 0.0) : Spectrum(1.0);
			std::vector<std::pair<Vec2, Vec2>> segments = geometry->segments();
			for (const std::pair<Vec2, Vec2> & segment : segments)
			{
				this->drawLine(Ivec2(rasterFromWorld2(segment.first)), Ivec2(rasterFromWorld2(segment.second)), false, color);
			}
		}
	}

	Vec2 rasterFromWorld2(const Vec2 & u)
	{
		return (u - mBound2.pMin) / (mBound2.pMax - mBound2.pMin) * Vec2(mFilm->mResolution);
	}

	Bound2 mBound2;
	std::atomic<Uint> mNumDrawnPaths;
	shared_ptr<Film> mFilm;
	shared_ptr<const Camera2> mCamera2;
	shared_ptr<const Scene> mScene;
};
