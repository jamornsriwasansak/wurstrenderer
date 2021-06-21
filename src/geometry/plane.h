#pragma once

#include "common/wurst.h"

#include "common/geometry.h"

struct PlaneGeometry : public Geometry
{
	PlaneGeometry(const Vec3 & origin,
				  const Vec3 & edge1,
				  const Vec3 & edge2,
				  const shared_ptr<Light> & areaLight):
		mOrigin(origin),
		mEdge1(edge1),
		mEdge2(edge2),
		Geometry(nullptr, areaLight)
	{
		mNormal = Math::Normalize(Math::Cross(mEdge1, mEdge2));
	}

	PlaneGeometry(const Vec3 & origin,
				  const Vec3 & edge1,
				  const Vec3 & edge2):
		PlaneGeometry(origin, edge1, edge2, nullptr)
	{}

	double area() const override
	{
		return Math::Length(mEdge1) * Math::Length(mEdge2);
	}

	Bound3 bound() const override
	{
		Bound3 result;

		const Vec3 a = mOrigin - mEdge1 * 0.5 - mEdge2 * 0.5;
		const Vec3 b = mOrigin + mEdge1 * 0.5 - mEdge2 * 0.5;
		const Vec3 c = mOrigin - mEdge1 * 0.5 + mEdge2 * 0.5;
		const Vec3 d = mOrigin + mEdge1 * 0.5 + mEdge2 * 0.5;

		result = Bound3::Union(result, a);
		result = Bound3::Union(result, b);
		result = Bound3::Union(result, c);
		result = Bound3::Union(result, d);

		// pad a little bit prevent bbox having 0 volume
		result.pMin -= Vec3(Math::SmallValue);
		result.pMax += Vec3(Math::SmallValue);
		
		return result;
	}

	Bound3 bound(const int primitiveIndex) const override
	{
		assert(primitiveIndex == 0);
		return bound();
	}

	double evalPdfA() const override
	{
		return 1.0 / area();
	}

	int getNumPrimitives() const override
	{
		return 1;
	}

	bool intersect(SimpleVertex * isect, Ray3 * ray, const int primitiveIndex) const override
	{
		assert(primitiveIndex == 0);

		const Vec3 p = Math::Cross(ray->direction, mEdge2);
		const double det = Math::Dot(mEdge1, p);
		if (det > -Math::SmallValue && det < Math::SmallValue) return false;

		const double invDet = 1.0 / det;
		const Vec3 t = ray->origin - mOrigin;

		const double u = Math::Dot(t, p) * invDet;
		if (u < -0.5 || u > 0.5) return false;

		const Vec3 q = Math::Cross(t, mEdge1);
		const double v = Math::Dot(ray->direction, q) * invDet;
		if (v < -0.5 || v > 0.5) return false;

		const double rayt = Math::Dot(mEdge2, q) * invDet;

		if (rayt < ray->tMin || rayt > ray->tMax) return false;

		ray->tMax = rayt;
		isect->mPosition = mOrigin + u * mEdge1 + v * mEdge2;
		isect->mGeometryNormal = mNormal;
		isect->mTextureCoord = Vec2(u, v) + Vec2(0.5);

		return true;
	}

	void precomputeForSampling() override
	{
	}

	SimpleVertex sample(double * pdfA, const Vec2 & sample) const override
	{
		if (pdfA) *pdfA = 1.0 / Math::ParallelogramArea(mEdge1, mEdge2);

		const Vec2 centeredSample = (sample - Vec2(0.5));

		SimpleVertex result;
		result.mGeometryNormal = mNormal;
		result.mPosition = mOrigin + centeredSample[0] * mEdge1 + centeredSample[1] * mEdge2;
		result.mTextureCoord = sample;

		return result;
	}

	void transform(const Matrix4 & matrix) override
	{
		mOrigin = matrix.transformPoint(mOrigin);
		mEdge1 = matrix.transformDirection(mEdge1);
		mEdge2 = matrix.transformDirection(mEdge2);
		mNormal = Math::Normalize(Math::Cross(mEdge1, mEdge2));
	}

	Vec3 mOrigin;

	// these three are bases of the plane. these bases aren't needed to be normalized.
	Vec3 mEdge1, mEdge2;
	// normalize(mEdge1 x mEdge2) = mNormal
	Vec3 mNormal;
};
