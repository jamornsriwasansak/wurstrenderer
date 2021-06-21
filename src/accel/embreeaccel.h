#pragma once

#include <map>
#include "common/wurst.h"
#include "common/path.h"

#include "embree3/rtcore.h"

#include "geometry/cube.h"
#include "geometry/plane.h"
#include "geometry/phongmesh.h"
#include "geometry/trianglemesh.h"

// TODO: implement the base class of this
struct EmbreeAccel
{
	struct EmbreeVertex { float x, y, z; };

	struct EmbreeTriangle { unsigned int v0, v1, v2; };

	struct EmbreeGeometryWrapper
	{
		EmbreeGeometryWrapper(shared_ptr<const Geometry> pg, const unsigned int geomId, const RTCGeometryType type):
			pg(pg), geomId(geomId), geomType(type) {}

		shared_ptr<const Geometry> pg;
		unsigned int geomId;
		RTCGeometryType geomType;
	};

	static void EmbreeErrorFunc(void* userPtr, RTCError code, const char* str)
	{
		switch (code)
		{
		case RTC_ERROR_NONE:
			std::cout << "RTC_ERROR_NONE" << std::endl;
			break;
		case RTC_ERROR_UNKNOWN:
			std::cout << "RTC_ERROR_UNKNOWN" << std::endl;
			break;
		case RTC_ERROR_INVALID_ARGUMENT:
			std::cout << "RTC_ERROR_INVALID_ARGUMENT" << std::endl;
			break;
		case RTC_ERROR_INVALID_OPERATION:
			std::cout << "RTC_ERROR_INVALID_OPERATION" << std::endl;
			break;
		case RTC_ERROR_OUT_OF_MEMORY:
			std::cout << "RTC_ERROR_OUT_OF_MEMORY" << std::endl;
			break;
		case RTC_ERROR_UNSUPPORTED_CPU:
			std::cout << "RTC_ERROR_UNSUPPORTED_CPU" << std::endl;
			break;
		case RTC_ERROR_CANCELLED:
			std::cout << "RTC_ERROR_CANCELLED" << std::endl;
			break;
		default:
			break;
		}
		std::cout << "EmbreeStr: " << str << std::endl;
		assert(false);
	}

	static void EmbreeUserBoundFunc(const struct RTCBoundsFunctionArguments * args)
	{
		const EmbreeGeometryWrapper * geom = reinterpret_cast<const EmbreeGeometryWrapper*>(args->geometryUserPtr);
		RTCBounds * bounds_o = args->bounds_o;
		int primId = int(args->primID);
		Bound3 b = geom->pg->bound(primId);
		bounds_o->lower_x = float(b.pMin[0]);
		bounds_o->lower_y = float(b.pMin[1]);
		bounds_o->lower_z = float(b.pMin[2]);
		bounds_o->upper_x = float(b.pMax[0]);
		bounds_o->upper_y = float(b.pMax[1]);
		bounds_o->upper_z = float(b.pMax[2]);
	}

	static void EmbreeUserIntersectFunc(const RTCIntersectFunctionNArguments* args)
	{
		const int* valid = args->valid;
		void* ptr = args->geometryUserPtr;
		RTCRayHit * rtcRayHit = reinterpret_cast<RTCRayHit*>(args->rayhit);
		int primId = int(args->primID);

		assert(args->N == 1);
		if (!valid[0]) return;

		Ray3 ray;
		ray.origin[0] = rtcRayHit->ray.org_x;
		ray.origin[1] = rtcRayHit->ray.org_y;
		ray.origin[2] = rtcRayHit->ray.org_z;
		ray.direction[0] = rtcRayHit->ray.dir_x;
		ray.direction[1] = rtcRayHit->ray.dir_y;
		ray.direction[2] = rtcRayHit->ray.dir_z;
		ray.tMin = rtcRayHit->ray.tnear;
		ray.tMax = rtcRayHit->ray.tfar;

		const EmbreeGeometryWrapper * geomWrapper = reinterpret_cast<const EmbreeGeometryWrapper*>(args->geometryUserPtr);
		SimpleVertex v;
		if (!geomWrapper->pg->intersect(&v, &ray, primId)) return;

		RTCHit potentialHit;
		potentialHit.u = float(v.mTextureCoord[0]);
		potentialHit.v = float(v.mTextureCoord[1]);
		potentialHit.instID[0] = args->context->instID[0];
		potentialHit.primID = args->primID;
		potentialHit.geomID = (unsigned int)(geomWrapper->geomId);
		potentialHit.Ng_x = float(v.mGeometryNormal[0]);
		potentialHit.Ng_y = float(v.mGeometryNormal[1]);
		potentialHit.Ng_z = float(v.mGeometryNormal[2]);

		RTCFilterFunctionNArguments fargs;
		int imask = -1;
		fargs.valid = &imask;
		fargs.geometryUserPtr = ptr;
		fargs.context = args->context;
		fargs.ray = reinterpret_cast<RTCRayN*>(args->rayhit);
		fargs.hit = reinterpret_cast<RTCHitN*>(&potentialHit);
		fargs.N = 1;

		rtcFilterIntersection(args, &fargs);

		// TODO:: alpha test
		//if (imask == -1)
		{
			rtcRayHit->ray.tfar = (float)ray.tMax;
			rtcRayHit->hit = potentialHit;
		}
	}

	static void EmbreeUserOcclusionFunc(const RTCOccludedFunctionNArguments* args)
	{
		const int* valid = args->valid;
		void* ptr = args->geometryUserPtr;
		RTCRay * rtcRay = reinterpret_cast<RTCRay*>(args->ray);
		int primId = int(args->primID);

		assert(args->N == 1);
		if (!valid[0]) return;

		Ray3 ray;
		ray.origin[0] = rtcRay->org_x;
		ray.origin[1] = rtcRay->org_y;
		ray.origin[2] = rtcRay->org_z;
		ray.direction[0] = rtcRay->dir_x;
		ray.direction[1] = rtcRay->dir_y;
		ray.direction[2] = rtcRay->dir_z;
		ray.tMin = rtcRay->tnear;
		ray.tMax = rtcRay->tfar;

		const EmbreeGeometryWrapper * geomWrapper = (const EmbreeGeometryWrapper *)args->geometryUserPtr;
		SimpleVertex v;
		if (!geomWrapper->pg->intersect(&v, &ray, primId)) return;

		RTCHit potentialHit;
		potentialHit.u = float(v.mTextureCoord[0]);
		potentialHit.v = float(v.mTextureCoord[1]);
		potentialHit.instID[0] = args->context->instID[0];
		potentialHit.primID = args->primID;
		potentialHit.geomID = (unsigned int)(geomWrapper->geomId);
		potentialHit.Ng_x = float(v.mGeometryNormal[0]);
		potentialHit.Ng_y = float(v.mGeometryNormal[1]);
		potentialHit.Ng_z = float(v.mGeometryNormal[2]);

		RTCFilterFunctionNArguments fargs;
		int imask = -1;
		fargs.valid = &imask;
		fargs.geometryUserPtr = ptr;
		fargs.context = args->context;
		fargs.ray = args->ray;
		fargs.hit = (RTCHitN*)&potentialHit;
		fargs.N = 1;

		rtcFilterOcclusion(args, &fargs);

		// TODO:: alpha test
		//if (imask == -1)
		{
			rtcRay->tfar = -std::numeric_limits<float>::infinity();
		}
	}

	static void EmbreeTriangleIntersectionFilterFunc(const RTCFilterFunctionNArguments * args)
	{
		assert(args->N == 1);
		if (args->context == nullptr) return;
		int * valid = args->valid;

		// ignore inactive rays 
		if (valid[0] != -1) return;

		RTCRayN * ray = args->ray;
		const TriangleMeshGeometry * triangleMeshGeometry = reinterpret_cast<const TriangleMeshGeometry *>(args->geometryUserPtr);

		// ignore hit if completely transparent 
		RTCHit * hit = reinterpret_cast<RTCHit*>(args->hit);
		if (triangleMeshGeometry->transparency(hit->primID, Vec2(hit->u, hit->v)) >= 0.5)
		{
			valid[0] = 0;
		}
	}

	static void EmbreeTriangleOcclusionFilterFunc(const RTCFilterFunctionNArguments * args)
	{
		assert(args->N == 1);
		if (args->context == nullptr) return;
		int * valid = args->valid;

		// ignore inactive rays 
		if (valid[0] != -1) return;

		RTCRayN * ray = args->ray;
		const TriangleMeshGeometry * triangleMeshGeometry = reinterpret_cast<const TriangleMeshGeometry *>(args->geometryUserPtr);

		// ignore hit if completely transparent 
		RTCHit * hit = reinterpret_cast<RTCHit*>(args->hit);
		if (triangleMeshGeometry->transparency(hit->primID, Vec2(hit->u, hit->v)) >= 0.5)
		{
			valid[0] = 0;
		}
	}

	EmbreeAccel(std::vector<shared_ptr<const TriangleMeshGeometry>> & geometries):
		EmbreeAccel()
	{
		for (shared_ptr<const TriangleMeshGeometry> & geometry : geometries) { addGeometry(geometry); }
		commit();
	}

	EmbreeAccel()
	{
		mRtcDevice = rtcNewDevice(nullptr);
		assert(mRtcDevice != nullptr);

		mRtcScene = rtcNewScene(mRtcDevice);
		assert(mRtcScene != nullptr);

		rtcSetDeviceErrorFunction(mRtcDevice, EmbreeErrorFunc, nullptr);

		rtcInitIntersectContext(&mIntersectContext);
	}

	void addUserGeometry(const shared_ptr<const Geometry> & p)
	{
		RTCGeometry geom = rtcNewGeometry(mRtcDevice, RTC_GEOMETRY_TYPE_USER);
		unsigned int geomId = rtcAttachGeometry(mRtcScene, geom);
		shared_ptr<EmbreeGeometryWrapper> egw = make_shared<EmbreeGeometryWrapper>(p, geomId, RTC_GEOMETRY_TYPE_USER);

		rtcSetGeometryUserPrimitiveCount(geom, p->getNumPrimitives());
		rtcSetGeometryUserData(geom, reinterpret_cast<void*>(egw.get()));
		rtcSetGeometryBoundsFunction(geom, EmbreeUserBoundFunc, nullptr);
		rtcSetGeometryIntersectFunction(geom, EmbreeUserIntersectFunc);
		rtcSetGeometryOccludedFunction(geom, EmbreeUserOcclusionFunc);

		rtcCommitGeometry(geom);
		rtcReleaseGeometry(geom);
		mIndexToGeometry[geomId] = egw;
	}

	void addGeometry(const shared_ptr<const CubeGeometry> & p)
	{
		addUserGeometry(p);
	}

	void addGeometry(const shared_ptr<const PlaneGeometry> & p)
	{
		addUserGeometry(p);
	}

	void addGeometry(const shared_ptr<const PhongMeshGeometry> & p)
	{
		addUserGeometry(p);
	}

	void addGeometry(const shared_ptr<const TriangleMeshGeometry> & p)
	{
		if (p->mTriangles.size() == 0) return;

		RTCGeometry geom = rtcNewGeometry(mRtcDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
		EmbreeVertex * vertices = reinterpret_cast<EmbreeVertex*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(EmbreeVertex), p->mPositions.size()));
		assert(vertices != nullptr);
		for (int i = 0; i < p->mPositions.size(); i++)
		{
			vertices[i].x = float(p->mPositions[i][0]);
			vertices[i].y = float(p->mPositions[i][1]);
			vertices[i].z = float(p->mPositions[i][2]);
		}

		EmbreeTriangle * triangles = reinterpret_cast<EmbreeTriangle*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(EmbreeTriangle), p->mTriangles.size()));
		assert(triangles != nullptr);
		for (int i = 0; i < p->mTriangles.size(); i++)
		{
			triangles[i].v0 = (unsigned int)(p->mTriangles[i][0]);
			triangles[i].v1 = (unsigned int)(p->mTriangles[i][1]);
			triangles[i].v2 = (unsigned int)(p->mTriangles[i][2]);
		}

		rtcSetGeometryVertexAttributeCount(geom, 1);
		float * textureCoords = reinterpret_cast<float*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT2, sizeof(float) * 2, p->mTextureCoords.size()));
		assert(textureCoords != nullptr);
		for (int i = 0; i < p->mTextureCoords.size(); i++)
		{
			textureCoords[i * 2 + 0] = float(p->mTextureCoords[i][0]);
			textureCoords[i * 2 + 1] = float(p->mTextureCoords[i][1]);
		}

		// allow alpha texture test if any
		rtcSetGeometryUserData(geom, const_cast<void*>(reinterpret_cast<const void*>(p.get())));
		rtcSetGeometryIntersectFilterFunction(geom, EmbreeTriangleIntersectionFilterFunc);
		rtcSetGeometryOccludedFilterFunction(geom, EmbreeTriangleOcclusionFilterFunc);

		rtcCommitGeometry(geom);
		unsigned int geomId = rtcAttachGeometry(mRtcScene, geom);
		mIndexToGeometry[geomId] = make_shared<EmbreeGeometryWrapper>(p, geomId, RTC_GEOMETRY_TYPE_TRIANGLE);
		rtcReleaseGeometry(geom);
	}

	Bound3 bound() const
	{
		RTCBounds bound;
		rtcGetSceneBounds(mRtcScene, &bound);
		Vec3 lower(bound.lower_x, bound.lower_y, bound.lower_z);
		Vec3 upper(bound.upper_x, bound.upper_y, bound.upper_z);
		return Bound3(lower, upper);
	}

	void boundingSphere(Vec3 * center, double * radius) const
	{
		RTCBounds bound;
		rtcGetSceneBounds(mRtcScene, &bound);
		Vec3 lower(bound.lower_x, bound.lower_y, bound.lower_z);
		Vec3 upper(bound.upper_x, bound.upper_y, bound.upper_z);
		if (center) *center = (lower + upper) * 0.5;
		if (radius) *radius = std::sqrt(Math::Distance2(lower, upper)) * 0.5;
	}

	bool isIntersect(const Ray3 & ray) const
	{
		RTCRay test;
		test.dir_x = float(ray.direction[0]);
		test.dir_y = float(ray.direction[1]);
		test.dir_z = float(ray.direction[2]);
		test.org_x = float(ray.origin[0]);
		test.org_y = float(ray.origin[1]);
		test.org_z = float(ray.origin[2]);
		test.tfar = float(ray.tMax);
		test.tnear = float(ray.tMin);
		test.time = 0.f;
		rtcOccluded1(mRtcScene, &mIntersectContext, &test);
		const bool hit = test.tfar < 0.0f;
		return hit;
	}

	shared_ptr<Vertex> intersect(Ray3 * ray) const
	{
		RTCRayHit query;
		query.ray.dir_x = float(ray->direction[0]);
		query.ray.dir_y = float(ray->direction[1]);
		query.ray.dir_z = float(ray->direction[2]);
		query.hit.geomID = RTC_INVALID_GEOMETRY_ID;
		query.ray.org_x = float(ray->origin[0]);
		query.ray.org_y = float(ray->origin[1]);
		query.ray.org_z = float(ray->origin[2]);
		query.hit.primID = RTC_INVALID_GEOMETRY_ID;
		query.ray.tfar = float(ray->tMax);
		query.ray.tnear = float(ray->tMin);
		query.ray.time = 0.f;

		rtcIntersect1(mRtcScene, &mIntersectContext, &query);
		const bool hit = (query.hit.geomID != RTC_INVALID_GEOMETRY_ID);
		if (!hit) return nullptr;

		// modify tmax
		ray->tMax = double(query.ray.tfar);

		// describe vertex information
		shared_ptr<const EmbreeGeometryWrapper> geometryWrapper = mIndexToGeometry.at(query.hit.geomID);
		shared_ptr<Vertex> vertex = nullptr;

		Vec2 textureCoord;
		if (geometryWrapper->geomType == RTC_GEOMETRY_TYPE_TRIANGLE)
		{
			float uv[2];
			rtcInterpolate0(rtcGetGeometry(mRtcScene, query.hit.geomID),
							query.hit.primID,
							query.hit.u,
							query.hit.v,
							RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE,
							0,
							uv,
							2);
			textureCoord = Vec2(uv[0], uv[1]);
		}
		else
		{
			textureCoord = Vec2(query.hit.u, query.hit.v);
		}

		if (geometryWrapper->pg->mAreaLight != nullptr)
		{
			shared_ptr<LightVertex> lv = make_shared<LightVertex>();
			lv->mGeometry = geometryWrapper->pg.get();
			lv->mLight = geometryWrapper->pg->mAreaLight.get();
			lv->mTextureCoord = textureCoord;

			vertex = lv;
		}
		else
		{
			shared_ptr<SurfaceVertex> sv = make_shared<SurfaceVertex>();
			sv->mBsdf = geometryWrapper->pg->mBsdf.get();
			sv->mGeometry = geometryWrapper->pg.get();
			sv->mTextureCoord = textureCoord;

			vertex = sv;
		}

		vertex->mPosition = ray->pMax();
		vertex->mGeometryNormal = Math::Normalize(Vec3(query.hit.Ng_x, query.hit.Ng_y, query.hit.Ng_z));
		vertex->mCoordFrame = CoordFrame3(vertex->mGeometryNormal);

		return vertex;
	}

	// implemented with purpose to imitate rasterization
	bool intersectTriangleInfo(TriangleInfoVertex * vertex, Ray3 * ray) const
	{
		RTCRayHit query;
		query.ray.dir_x = float(ray->direction[0]);
		query.ray.dir_y = float(ray->direction[1]);
		query.ray.dir_z = float(ray->direction[2]);
		query.hit.geomID = RTC_INVALID_GEOMETRY_ID;
		query.ray.org_x = float(ray->origin[0]);
		query.ray.org_y = float(ray->origin[1]);
		query.ray.org_z = float(ray->origin[2]);
		query.hit.primID = RTC_INVALID_GEOMETRY_ID;
		query.ray.tfar = float(ray->tMax);
		query.ray.tnear = float(ray->tMin);
		query.ray.time = 0.f;

		rtcIntersect1(mRtcScene, &mIntersectContext, &query);
		const bool hit = (query.hit.geomID != RTC_INVALID_GEOMETRY_ID);
		if (!hit) return false;

		ray->tMax = query.ray.tfar;

		auto mappedTriMesh = mIndexToGeometry.find(query.hit.geomID);
		if (mappedTriMesh == mIndexToGeometry.end()) return false;
		if (mappedTriMesh->second->geomType != RTC_GEOMETRY_TYPE_TRIANGLE) return false;
		const shared_ptr<const TriangleMeshGeometry> triMesh = static_pointer_cast<const TriangleMeshGeometry>(mappedTriMesh->second->pg);

		vertex->mBarycentricCoords = Vec2(query.hit.u, query.hit.v);
		vertex->mGeometryNormal = Math::Normalize(Vec3(query.hit.Ng_x, query.hit.Ng_y, query.hit.Ng_z));
		vertex->mPosition = ray->pMax();
		float uv[2];
		rtcInterpolate0(rtcGetGeometry(mRtcScene, query.hit.geomID), query.hit.primID, query.hit.u, query.hit.v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, uv, 2);
		vertex->mTextureCoords = Vec2(uv[0], uv[1]);
		vertex->mGeometry = triMesh;

		vertex->mTriangleId = query.hit.primID;
		vertex->mGeometryId = query.hit.geomID;
		for (int i = 0; i < 3; i++) { vertex->mIndices[i] = triMesh->mTriangles[vertex->mTriangleId][i]; }

		return true;
	}

	void commit()
	{
		rtcCommitScene(mRtcScene);
	}
	
	std::vector<shared_ptr<const EmbreeGeometryWrapper>> mUserGeomWrappers;
	std::map<unsigned int, shared_ptr<const EmbreeGeometryWrapper>> mIndexToGeometry;
	RTCDevice mRtcDevice = nullptr;
	RTCScene mRtcScene = nullptr;
	mutable RTCIntersectContext mIntersectContext;
};