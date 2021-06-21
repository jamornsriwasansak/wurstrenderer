#pragma once

#include "common/wurst.h"
#include "common/cdftable.h"
#include "common/geometry.h"
#include "common/texture.h"

struct TriangleMeshGeometry;

// for intersection test
struct TriangleInfoVertex
{
	Vec3 mPosition;
	Vec3 mGeometryNormal;
	Vec2 mTextureCoords;
	Vec2 mBarycentricCoords;
	shared_ptr<const TriangleMeshGeometry> mGeometry;

	int mTriangleId;
	int mGeometryId;
	int mIndices[3];
};

struct TriangleMeshGeometry : public Geometry
{
	static bool Intersect(double * b1Ptr, double * b2Ptr, double * tPtr, const Vec3 & p0, const Vec3 & p1, const Vec3 & p2, const Ray3 & ray)
	{
		double & u = *b1Ptr;
		double & v = *b2Ptr;
		double & t = *tPtr;

		Vec3 e1 = p1 - p0, e2 = p2 - p0;
		Vec3 pvec = Math::Cross(ray.direction, e2);
		double det = Math::Dot(e1, pvec);
		if (det == 0) return false;

		double invDet = double(1.0f / det);

		Vec3 tvec = ray.origin - p0;
		u = Math::Dot(tvec, pvec) * invDet;
		if (u < 0.0f || u > 1.0f) return false;

		Vec3 qvec = Math::Cross(tvec, e1);
		v = Math::Dot(ray.direction, qvec) * invDet;

		if (v < 0.0 || u + v > 1.0) return false; 

		t = Math::Dot(e2, qvec) * invDet;
		if (t < ray.tMin || t > ray.tMax) return false;

		return true;
	}

	TriangleMeshGeometry()
	{
	}

	TriangleMeshGeometry(const shared_ptr<Light> & areaLight) : Geometry(nullptr, areaLight)
	{
	}

	double area() const override
	{
		return mArea;
	}

	Bound3 bound() const override
	{
		Bound3 result;
		for (const Vec3 & position : mPositions)
		{
			result = Bound3::Union(result, position);
		}
		return result;
	}

	Bound3 bound(const int primitiveIndex) const override
	{
		Bound3 result;
		Ivec3 indices = mTriangles[primitiveIndex];
		for (int i = 0; i < 3; i++)
		{
			result = Bound3::Union(result, mPositions[indices[i]]);
		}
		return result;
	}

	double evalPdfA() const override
	{
		return mInvArea;
	}

	int getNumPrimitives() const override
	{
		return int(mTriangles.size());
	}

	bool intersect(SimpleVertex * isect, Ray3 * ray, const int primitiveIndex) const override
	{
		double u, v, t;
		const Ivec3 tri = mTriangles[primitiveIndex];
		if (Intersect(&u, &v, &t, mPositions[tri[0]], mPositions[tri[1]], mPositions[tri[2]], *ray))
		{
			ray->tMax = t;
			return true;
		}
		return false;
	}

	void precomputeForSampling() override
	{
		// compute all triangle areas
		std::vector<int> triangleFaceIndex(mTriangles.size());
		std::vector<double> triangleAreas(mTriangles.size());
		double sumArea = 0.0;
		for (int iTriangle = 0; iTriangle < mTriangles.size(); iTriangle++)
		{
			triangleFaceIndex[iTriangle] = iTriangle;

			Vec3 positions[3];
			for (int iVertex = 0; iVertex < 3; iVertex++) { positions[iVertex] = mPositions[mTriangles[iTriangle][iVertex]]; }
			double area = Math::ComputeTriangleArea(positions[0], positions[1], positions[2]);
			sumArea += area;
			triangleAreas[iTriangle] = area;
		}
		mArea = sumArea;
		mInvArea = 1.0 / sumArea;
		mTriangleAreasTable = CdfTable(triangleAreas);
	}

	SimpleVertex sample(double * pdfA, const Vec2 & sample) const override
	{
		SimpleVertex result;

		// sample a triangle
		Vec2 lsample = sample;
		int triangleIndex = mTriangleAreasTable.sample(nullptr, &lsample[0], lsample[0]);

		// sample a barycentric coordinate
		Vec2 baryCoord = Mapping::BarycentricFromSquare(lsample);

		// map the barycentric coordinate to position on triangle
		Vec3 positions[3];
		Vec2 texCoords[3];
		for (int i = 0; i < 3; i++)
		{
			int vertexIndex = mTriangles[triangleIndex][i];
			positions[i] = mPositions[vertexIndex];
			texCoords[i] = mTextureCoords[vertexIndex];
		}

		result.mPosition = Math::BarycentricInterpolate(positions[0], positions[1], positions[2], baryCoord);
		result.mGeometryNormal = Math::Normalize(Math::Cross(positions[1] - positions[0], positions[2] - positions[0]));
		result.mTextureCoord = Math::BarycentricInterpolate(texCoords[0], texCoords[1], texCoords[2], baryCoord);

		if (pdfA) *pdfA = mInvArea;

		return result;
	}

	double transparency(const int triangleIndex, const Vec2 & uv) const override
	{
		if (mAlpha == nullptr) return 0.0;

		Ivec3 t = mTriangles[triangleIndex];
		const Vec2 texCoord = Math::BarycentricInterpolate(mTextureCoords[t[0]],
														   mTextureCoords[t[1]],
														   mTextureCoords[t[2]],
														   uv);
		return 1.0 - mAlpha->eval(texCoord);
	}

	void transform(const Matrix4 & matrix) override
	{
		for (int i = 0; i < mPositions.size(); i++)
		{
			mPositions[i] = matrix.transformPoint(mPositions[i]);
		}
	}

	// 2d

	Bound2 bound2() const override
	{
		Bound2 result;
		for (const Vec2 & position : mPositions2)
		{
			result = Bound2::Union(result, position);
		}
		return result;
	}

	Bound2 bound2(const int primitiveIndex) const override
	{
		Bound2 result;
		Ivec3 indices = mTriangles[primitiveIndex];
		for (int i = 0; i < 3; i++)
		{
			result = Bound2::Union(result, mPositions2[indices[i]]);
		}
		return result;
	}

	void precomputeForSampling2(const double planeZ) override
	{
		mPositions2 = std::vector<Vec2>();
		mTextureCoords2 = std::vector<Vec2>();
		mSegments = std::vector<Ivec2>(mTriangles.size(), Ivec2(0));
		double sumLength = 0.0;
		std::vector<double> lineLengths(mTriangles.size(), 0.0);

		// extract lines from intersection between planez and triangles
		for (int iTriangle = 0; iTriangle < mTriangles.size(); iTriangle++)
		{
			Vec3 positions[3];
			Vec2 texCoords[3];
			for (int i = 0; i < 3; i++)
			{
				int vertexIndex = mTriangles[iTriangle][i];
				positions[i] = mPositions[vertexIndex];
				texCoords[i] = mTextureCoords[vertexIndex];
			}

			std::vector<Vec2> barycentricCoords;
			std::vector<Vec3> segmentPositions;
			std::vector<Vec2> segmentTextureCoords;
			for (int i = 0; i < 3; i++)
			{
				double t;
				if (Math::IntersectPlaneSegment(&t, planeZ, positions[i], positions[(i + 1) % 3]))
				{
					Vec3 barycentricCoord(0.0);
					barycentricCoord[i] = t;
					barycentricCoord[(i + 1) % 3] = 1.0 - t;
					barycentricCoords.push_back(Vec2(barycentricCoord[0], barycentricCoord[1]));

					Vec3 segmentPosition = (positions[(i + 1) % 3] - positions[i]) * t + positions[i];
					segmentPositions.push_back(segmentPosition);

					Vec2 segmentTexCoords = (texCoords[(i + 1) % 3] - texCoords[i]) * t + texCoords[i];
					segmentTextureCoords.push_back(segmentTexCoords);
				}
			}

			if (barycentricCoords.size() == 2)
			{
				int index = int(mPositions2.size());
				mPositions2.push_back(Vec2(segmentPositions[0][0], segmentPositions[0][1]));
				mPositions2.push_back(Vec2(segmentPositions[1][0], segmentPositions[1][1]));
				mTextureCoords2.push_back(segmentTextureCoords[0]);
				mTextureCoords2.push_back(segmentTextureCoords[1]);
				mSegments[iTriangle] = Ivec2(index, index + 1);

				const double length = Math::Distance(segmentPositions[0], segmentPositions[1]);
				sumLength += length;
				lineLengths[iTriangle] = length;
			}
		}

		mLength = sumLength;
		mInvLength = 1.0 / sumLength;
		mSegmentLengthsTable = CdfTable(lineLengths);
	}

	SimpleVertex2 sample2(double * pdfL, const double sample) const override
	{
		SimpleVertex2 result;

		// sample a triangle
		double lsample = sample;
		int segmentIndex = mSegmentLengthsTable.sample(nullptr, &lsample, lsample);

		// pull all the data
		Vec2 positions2[2];
		Vec2 texCoords2[2];
		for (int i = 0; i < 2; i++)
		{
			int vertexIndex = mSegments[segmentIndex][i];
			positions2[i] = mPositions2[vertexIndex];
			texCoords2[i] = mTextureCoords2[vertexIndex];
		}

		// compute geometry normal
		Vec2 geomNormal = positions2[0] - positions2[1];
		geomNormal = Math::Normalize(Vec2(geomNormal[1], geomNormal[0]));

		// record
		result.mPosition2 = Math::Interpolate(positions2[0], positions2[1], sample);
		result.mGeometryNormal2 = geomNormal;
		result.mTextureCoord = Math::Interpolate(texCoords2[0], texCoords2[1], sample);

		if (pdfL) *pdfL = mInvLength;

		return result;
	}

	std::vector<std::pair<Vec2, Vec2>> segments() const override
	{
		std::vector<std::pair<Vec2, Vec2>> result;
		for (int iSegment = 0; iSegment < mSegments.size(); iSegment++)
		{
			int iPos0 = mSegments[iSegment][0];
			int iPos1 = mSegments[iSegment][1];

			// the segement is degenerate
			if (iPos0 == iPos1) continue;

			result.push_back(make_pair(mPositions2[iPos0], mPositions2[iPos1]));
		}
		return result;
	}

	// arguments for 3d

	std::vector<Vec3> mPositions;
	std::vector<Vec3> mNormals;
	std::vector<Vec2> mTextureCoords;
	std::vector<Ivec3> mTriangles;
	double mInvArea = 0.0;
	double mArea = 0.0;
	CdfTable mTriangleAreasTable;

	// arguments for 2d

	std::vector<Vec2> mPositions2;
	std::vector<Vec2> mTextureCoords2;
	std::vector<Ivec2> mSegments;
	double mLength = 0.0;
	double mInvLength = 0.0;
	CdfTable mSegmentLengthsTable;
};

