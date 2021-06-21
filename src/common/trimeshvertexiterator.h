#pragma once

#include "common/scene.h"
#include "common/parallel.h"
#include "common/path.h"

// this struct is to replicate fixed pipeline of OpenGL / DirectX
struct TriangleMeshVerticesVisitor
{
	TriangleMeshVerticesVisitor(const shared_ptr<const Scene> & scene) :
		mScene(scene)
	{
		// check if all meshes are triangle mesh
		for (int iMesh = 0; iMesh < mScene->mGeometries.size(); iMesh++)
		{
			if (shared_ptr<const TriangleMeshGeometry> triMesh = dynamic_pointer_cast<const TriangleMeshGeometry>(mScene->mGeometries[iMesh])) {} // do nothing
			else { throw std::runtime_error("prt doesn't accept non trianglemeshes"); }
		}

		// init num points
		mNumPoints = std::vector<int>(mScene->mGeometries.size());
		forEachMesh([&](int iMesh, shared_ptr<const TriangleMeshGeometry> & triMesh)
		{
			mNumPoints[iMesh] = int(triMesh->mPositions.size());
		});

		// init cumulative num points
		mCumulativeNumPoints = std::vector<int>(mScene->mGeometries.size());
		mCumulativeNumPoints[0] = 0;
		forEachMesh([&](int iMesh, shared_ptr<const TriangleMeshGeometry> & triMesh)
		{
			if (iMesh == 0) return;
			mCumulativeNumPoints[iMesh] += mNumPoints[iMesh - 1] + mCumulativeNumPoints[iMesh - 1];
		});
		mNumTotalPoints = mCumulativeNumPoints.back() + mNumPoints[mScene->mGeometries.size() - 1];
	}

	void forEachMesh(const std::function<void(int iMesh, shared_ptr<const TriangleMeshGeometry> & triMesh)> & callbackFunc) const
	{
		for (int iMesh = 0; iMesh < mScene->mGeometries.size(); iMesh++)
			if (shared_ptr<const TriangleMeshGeometry> triMesh = dynamic_pointer_cast<const TriangleMeshGeometry>(mScene->mGeometries[iMesh]))
				callbackFunc(iMesh, triMesh);
	}

	template <typename ParallelMode=Serial>
	void forEachPoint(const std::function<void(int iMesh, shared_ptr<const TriangleMeshGeometry> & triMesh, int iPoint)> & callbackFunc) const
	{
		forEachPoint<ParallelMode>(0, mNumTotalPoints, callbackFunc);
	}

	template <typename ParallelMode=Serial>
	void forEachPoint(const int start, const int end, const std::function<void(int iMesh, shared_ptr<const TriangleMeshGeometry> & triMesh, int iPoint)> & callbackFunc) const
	{
		ParallelMode::Split(start, end, [&] (const int iStart, const int iEnd) 
		{
			for (int gPoint = iStart; gPoint < iEnd; gPoint++)
			{
				int iMesh, iPoint;
				getIMeshIPoint(&iMesh, &iPoint, gPoint);
				shared_ptr<const TriangleMeshGeometry> triMesh = static_pointer_cast<const TriangleMeshGeometry>(mScene->mGeometries[iMesh]);
				callbackFunc(iMesh, triMesh, iPoint);
			}
		});
	}

	void getIMeshIPoint(int * iMesh, int * iPoint, const int index) const
	{
		// find iPoint and iMesh
		auto iter = std::upper_bound(mCumulativeNumPoints.begin(), mCumulativeNumPoints.end(), index);
		--iter;
		if (iMesh) *iMesh = int(std::distance(mCumulativeNumPoints.begin(), iter));
		if (iPoint) *iPoint = index - *iter;
	}

	int index(int iMesh, int iPoint) const
	{
		return mCumulativeNumPoints[iMesh] + iPoint;
	}

	shared_ptr<SurfaceVertex> getSurfaceVertex(const int iMesh, const int iPoint) const
	{
		shared_ptr<const TriangleMeshGeometry> triMesh = static_pointer_cast<const TriangleMeshGeometry>(mScene->mGeometries[iMesh]);

		const Vec3 position = triMesh->mPositions[iPoint];
		const Vec3 normal = triMesh->mNormals[iPoint];
		const Vec2 texCoord = triMesh->mTextureCoords[iPoint];

		// initialize surface vertex
		shared_ptr<SurfaceVertex> sv = make_shared<SurfaceVertex>();
		sv->mBsdf = triMesh->mBsdf.get();
		sv->mCoordFrame = CoordFrame3(normal);
		sv->mGeometry = triMesh.get();
		sv->mGeometryNormal = normal;
		sv->mPathContrib = Spectrum(1.0);
		sv->mPosition = position;
		sv->mSubpathSampler = nullptr;
		sv->mTextureCoord = texCoord;

		return sv;
	}

	shared_ptr<SurfaceVertex> getSurfaceVertex(const int globalPointIndex) const
	{
		int iMesh, iPoint;
		getIMeshIPoint(&iMesh, &iPoint, globalPointIndex);
		return getSurfaceVertex(iMesh, iPoint);
	}
	
	shared_ptr<const Scene> mScene;
	std::vector<int> mNumPoints;
	std::vector<int> mCumulativeNumPoints;
	int mNumTotalPoints;
};
