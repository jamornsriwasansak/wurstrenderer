#pragma once

#include "common/wurst.h"

struct SimpleVertex
{
	Vec3 mPosition;
	Vec3 mGeometryNormal;
	Vec2 mTextureCoord;
};

struct SimpleVertex2
{
	Vec2 mPosition2;
	Vec2 mGeometryNormal2;
	Vec2 mTextureCoord;
};

struct Geometry
{
	Geometry()
	{
	}

	Geometry(const shared_ptr<const Bsdf> & bsdf, const shared_ptr<const Light> & areaLight) :
		mBsdf(bsdf), mAreaLight(areaLight)
	{
	}

	virtual double area() const = 0;

	virtual Bound3 bound() const = 0;

	virtual Bound3 bound(const int primitiveIndex) const = 0;

	virtual double evalPdfA() const = 0;

	virtual int getNumPrimitives() const = 0;

	virtual bool intersect(SimpleVertex * isect, Ray3 * ray, const int primitiveIndex) const = 0;

	virtual void precomputeForSampling() = 0;

	virtual SimpleVertex sample(double * pdfA, const Vec2 & sample) const = 0;

	virtual double transparency(const int triangleIndex, const Vec2 & uv) const { throw std::runtime_error("unimpl"); }

	virtual void transform(const Matrix4 & matrix) = 0;

	// 2d stuffs

	virtual Bound2 bound2() const { throw std::runtime_error("unimpl"); }

	virtual Bound2 bound2(const int primitiveIndex) const { throw std::runtime_error("unimpl"); }

	virtual void precomputeForSampling2(const double planeZ) { throw std::runtime_error("unimpl"); };

	virtual SimpleVertex2 sample2(double * pdfL, const double sample) const { throw std::runtime_error("unimpl"); };

	virtual std::vector<std::pair<Vec2, Vec2>> segments() const { throw std::runtime_error("unimpl"); }

	shared_ptr<const Light> mAreaLight = nullptr;
	shared_ptr<const Texture<double>> mAlpha = nullptr;
	shared_ptr<const Bsdf> mBsdf = nullptr;

	shared_ptr<const Medium> mMediumFrontface = nullptr;
	shared_ptr<const Medium> mMediumBackface = nullptr;
};