#pragma once

#include "common/wurst.h"

#include "common/film.h"

struct Camera2
{
	Camera2(const Ivec2 & resolution): mFilm(make_shared<Film>(resolution)) {}

	virtual Bound2 bound2() const = 0;

	virtual Spectrum samplePerPixelWe0ByPdf(CameraVertex * cv, double * pdfL, const Ivec2 & pixel, const double sample) const = 0;

	virtual Spectrum sampleWe0ByPdf(CameraVertex * cv, Vec2 * raster, double * pdfL, const double sample) const { throw std::runtime_error("unimpl"); }

	virtual Spectrum sampleWe1CosByPdf(Vec2 * direction, double * pdfW, const CameraVertex & cv, const double sample) const = 0;

	virtual std::vector<std::pair<Vec2, Vec2>> segments() const = 0;

	shared_ptr<Film> mFilm = nullptr;
};
