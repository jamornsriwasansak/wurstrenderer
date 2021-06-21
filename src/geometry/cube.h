#pragma once

#include "common/wurst.h"

#include "common/geometry.h"

struct CubeGeometry : public Geometry
{
	CubeGeometry()
	{}

	double area() const override
	{
		return (mArea[0] + mArea[1] + mArea[2]) * 2;
	}

	Bound3 bound() const override
	{
		Bound3 result;
		result = Bound3::Union(result, mToGlobal.transformPoint(Vec3(0.5, 0.5, 0.5)));
		result = Bound3::Union(result, mToGlobal.transformPoint(Vec3(-0.5, 0.5, 0.5)));
		result = Bound3::Union(result, mToGlobal.transformPoint(Vec3(0.5, -0.5, 0.5)));
		result = Bound3::Union(result, mToGlobal.transformPoint(Vec3(0.5, 0.5, -0.5)));
		result = Bound3::Union(result, mToGlobal.transformPoint(Vec3(0.5, -0.5, -0.5)));
		result = Bound3::Union(result, mToGlobal.transformPoint(Vec3(-0.5, 0.5, -0.5)));
		result = Bound3::Union(result, mToGlobal.transformPoint(Vec3(-0.5, -0.5, 0.5)));
		result = Bound3::Union(result, mToGlobal.transformPoint(Vec3(-0.5, -0.5, -0.5)));
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

		// transform ray 
		Vec3 rayOrigin = mToLocal.transformPoint(ray->origin);
		Vec3 rayDirection = mToLocal.transformDirection(ray->direction);

		// intersect unit cube
		double t0 = ray->tMin;
		double t1 = ray->tMax;
		for (int iAxis = 0; iAxis < 3; iAxis++)
		{
			if (rayDirection[iAxis] == 0)
			{
				if (rayOrigin[iAxis] < -0.5 || rayOrigin[iAxis] > 0.5) { return false; }
			}
			else
			{
				double invRayDir = 1.f / rayDirection[iAxis];
				double tNear = (-0.5 - rayOrigin[iAxis]) * invRayDir;
				double tFar = (0.5 - rayOrigin[iAxis]) * invRayDir;
				if (tNear > tFar) { std::swap(tNear, tFar); }

				t0 = std::max(tNear, t0);
				t1 = std::min(tFar, t1);
				if (t0 > t1) return false;
			}
		}

		Vec3 position;
		Vec3 normal(0.0);
		Vec2 texCoord;

		if (t0 > ray->tMin && t0 < ray->tMax) // intersect from outside the box
		{
			// position
			position = rayOrigin + t0 * rayDirection;

			// normal
			int maxDim = int(Math::MaxDimension(Math::Abs(position)));
			normal[maxDim] = position[maxDim] > 0.0 ? 1.0 : -1.0;

			// tex coord
			texCoord = Vec2(0.0);

			ray->tMax = t0;
		}
		else // intersect from inside
		{
			// position
			position = rayOrigin + t1 * rayDirection;

			// normal
			int maxDim = int(Math::MaxDimension(Math::Abs(position)));
			normal[maxDim] = position[maxDim] > 0.0 ? -1.0 : 1.0;

			// tex coord
			texCoord = Vec2(0.0);

			ray->tMax = t1;
		}

		isect->mPosition = mToGlobal.transformPoint(position);
		isect->mGeometryNormal = mToGlobal.transformDirection(normal);
		isect->mTextureCoord = texCoord;
		return true;
	}

	void precomputeForSampling() override
	{
		// extract scaling matrix
		assert(Math::IsApprox(mToGlobal.cols[0][3], 0.0));
		assert(Math::IsApprox(mToGlobal.cols[1][3], 0.0));
		assert(Math::IsApprox(mToGlobal.cols[2][3], 0.0));

		mScale = Vec3(Math::Length(mToGlobal.cols[0]), Math::Length(mToGlobal.cols[1]), Math::Length(mToGlobal.cols[2]));

		mArea = std::vector<double>(6);
		for (int i = 0; i < 2; i++)
		{
			mArea[0 + i*3] = mScale[1] * mScale[2];
			mArea[1 + i*3] = mScale[0] * mScale[2];
			mArea[2 + i*3] = mScale[0] * mScale[1];
		}

		mAreaTable = CdfTable(mArea);
	}

	SimpleVertex sample(double * pdfA, const Vec2 & sample) const override
	{
		Vec2 lsample = sample;
		int dimFace = mAreaTable.sample(pdfA, &lsample[0], lsample[0]);
		double faceSign = (dimFace / 3 == 0) ? 1.0 : -1.0;

		SimpleVertex sv;

		// normal
		Vec3 normal(0.0);
		normal[dimFace % 3] = faceSign;
		sv.mGeometryNormal = mToGlobal.transformDirection(normal);

		// position
		Vec3 position(0.0);
		position[dimFace % 3] = faceSign * 0.5;
		position[(dimFace + 1) % 3] = lsample[0] - 0.5;
		position[(dimFace + 2) % 3] = lsample[1] - 0.5;
		sv.mPosition = mToGlobal.transformPoint(position);

		// texture coord
		sv.mTextureCoord = Vec2(0.0);

		return sv;
	}

	void transform(const Matrix4 & matrix) override
	{
		mToGlobal = matrix * mToGlobal ;
		mToLocal = Matrix4::Inverse(mToGlobal);
	}

	std::vector<double> mArea;
	CdfTable mAreaTable;
	Vec3 mScale; // scale component extracted from mToGlobal Matrix
	Matrix4 mToLocal = Matrix4();
	Matrix4 mToGlobal = Matrix4();
};
