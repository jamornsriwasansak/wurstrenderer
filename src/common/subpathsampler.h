#pragma once

#include "common/wurst.h"
#include "common/util/file.h"

struct Subpath
{
	std::vector<shared_ptr<Vertex>> mVertices;
	std::vector<Spectrum> mContribs;
	std::vector<double> mPdfAs;
};

struct SubpathSampler
{
	enum class Request : uint8_t
	{
		None = 0,
		PdfA = 1 << 0
	};

	virtual double evalCameraPdfW(const CameraVertex & c, const Vec3 & w) const = 0;

	virtual double evalLightPdfA(const LightVertex & lv) const = 0;
};

DECLARE_ENUM_OPERATORS(SubpathSampler::Request);