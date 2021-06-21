#pragma once

#include "common/light.h"
#include "common/texture.h"
#include "common/sphericalrepresentation.h"
#include "common/path.h"
#include "texture/constanttexture.h"

struct EnvmapLight : public Light
{
	EnvmapLight(const shared_ptr<const SphRep<Spectrum>> & sphRep, const Vec3 & center, const double sceneRadius):
		mSphRep(sphRep)
	{
		updateSceneRadius(center, sceneRadius);
	}

	EnvmapLight(const shared_ptr<const SphRep<Spectrum>> & sphRep):
		EnvmapLight(sphRep, Vec3(0.0), 0.0)
	{
	}

	double cdfWeight() const override
	{
		// TODO:: this weight is under optimal
		return 1.23;
	}

	// return pdf in solid angle domain
	double evalPdfLe0(const LightVertex & vertex) const override
	{
        assert(Math::IsNormalized(vertex.mPosition));
		return mSphRep->pdfA(vertex.mPosition);
	}

	// return pdf in area domain
	double evalPdfLe1(const LightVertex & vertex, const Vec3 & outgoing) const override
	{
		return Math::Pi * mRadius * mRadius;
	}

	Spectrum evalLe0(const LightVertex & vertex) const override
	{
		assert(Math::IsNormalized(vertex.mPosition));
		return mSphRep->eval(vertex.mPosition);
	}

	Spectrum evalLe1(const LightVertex & vertex, const Vec3 & w) const override
	{
		return Spectrum(1.0);
	}

	shared_ptr<EnvmapVertex> intersect(Ray3 * ray) const
	{
		shared_ptr<EnvmapVertex> lv = make_shared<EnvmapVertex>();
		lv->mCoordFrame = CoordFrame3();
		lv->mGeometry = nullptr;
		lv->mLight = this;
		lv->mPosition = ray->direction;
		return lv;
	}
	
	// firstly sample the outgoing direction
	Spectrum sampleLe0ByPdf(LightVertex * lv, double * pdfW, const Vec2 & sample) const override
	{
		Vec3 position;
		double pdf;
		Spectrum result = mSphRep->sample(&position, &pdf, sample);

		if (lv)
		{
			lv->mCoordFrame = CoordFrame3();
			lv->mPosition = position;
			lv->mPathContrib = result;
			lv->mLight = this;
		}
		if (pdfW) *pdfW = pdf;

		return result;
	}

	// secondly sample the lightsource position
	Spectrum sampleLe1CosByPdf(Vec3 * outgoing, double * pdfA, const LightVertex & vertex, const Vec2 & sample) const override
	{
		assert(vertex.hasType(VertexType::Envmap));

		// sample a point on the disk of mRadius
		Vec2 p = Mapping::DiskFromSquare(sample) * mRadius;
		
		// translate the point along y  by -r
		Vec3 q(p[0], -mRadius, p[1]);

		// rotate the the point according to sampled direction
		const Vec3 & direction = -vertex.mPosition;
		const CoordFrame3 frame(direction);
		const double area = Math::Pi * Math::Square(mRadius);
		const double pdf = 1.0 / area;

		if (outgoing) *outgoing = frame.toWorld(q) + mCenter;
		if (pdfA) *pdfA = pdf;
		return Spectrum(area);
	}

	void updateSceneRadius(const Vec3 & center, const double sceneRadius)
	{
		mCenter = center;
		mRadius = sceneRadius;
	}

	shared_ptr<const SphRep<Spectrum>> mSphRep = nullptr;

	Vec3 mCenter;
	double mRadius = 0.0;
};