#pragma once

#include "common/camera2.h"
#include "common/scene.h"
#include "common/path.h"
#include "common/sampler.h"

#include "common/util/file.h"

struct SimpleSubpathSampler2
{
	SimpleSubpathSampler2(const shared_ptr<const Scene> & scene, const bool doStratification):
		mScene(scene),
		mDoStratification(doStratification)
	{
	}

	static void RandomWalk(Sampler * sampler,
						   const Scene & scene,
						   shared_ptr<Vertex> prevVertex,
						   Ray2 ray,
						   Spectrum prevContrib,
						   double prevPdfT,
						   const Uint numVertices,
						   const std::function<void(shared_ptr<Vertex> & currentVertex,
													shared_ptr<Vertex> prevVertex,
													const Spectrum & contrib,
													const double pdfW)> & callbackFunc)
	{
		assert(sampler != nullptr);

		if (numVertices == 0) return;

		Uint countVertices = 0;

		while (true)
		{
			Ray3 ray3(Vec3(ray.origin[0], ray.origin[1], scene.mPlaneZ), Vec3(ray.direction[0], ray.direction[1], 0.0));
			shared_ptr<Vertex> vertex = scene.intersect(&ray3);
			if (vertex == nullptr) return;

			// update dir to parent vertex, rectify geometry normal and basis
			vertex->mGeometryNormal2 = Math::Normalize(Vec2(vertex->mGeometryNormal[0], vertex->mGeometryNormal[1]));
			vertex->mCoordFrame2 = CoordFrame2(vertex->mGeometryNormal2);

			// callback 
			if (callbackFunc) callbackFunc(vertex, prevVertex, prevContrib, prevPdfT);

			// num vertices exceed, break!
			if (++countVertices == numVertices) break;

			if (vertex->getVertexType() == VertexType::Surface) // Hit surface
			{
				// init data for sampling next direction
				shared_ptr<SurfaceVertex> sv = Vertex::SafeCast<SurfaceVertex>(vertex);

				Vec2 inLocal = sv->mCoordFrame2.toLocal(-ray.direction);
				Vec2 outLocal;

				// sample next direction 
				Spectrum bsdfCosByPdf = sv->mBsdf->sampleBsdfCosByPdf2(&outLocal, nullptr, *sv, inLocal, sampler->get1d());
				if (Math::IsZero(bsdfCosByPdf)) { return; }
				assert(Math::IsNormalized(outLocal));
				prevContrib = bsdfCosByPdf;

				// update ray
				ray = Ray2(sv->mPosition2, sv->mCoordFrame2.toWorld(outLocal));
				prevVertex = vertex;
			}
			else if (vertex->getVertexType() == VertexType::Medium)
			{
				// TODO:: implement
			}
			else // Hit light or camera
			{
				return;
			}
		}
	}

	void createCameraSubpath(Sampler * sampler,
							 Vec2 * raster,
							 const Camera2 & camera,
							 const Ivec2 & pixel,
							 const Uint numVertices,
							 const std::function<void(shared_ptr<CameraVertex> endpointVertex,
													  const Spectrum & contrib,
													  const double pdfL)> & epCallbackFunc,
							 const std::function<void(shared_ptr<Vertex> currentVertex,
													  shared_ptr<Vertex> prevVertex,
													  const Spectrum & contrib,
													  const double pdfL)> & callbackFunc = nullptr) const
	{
		assert(sampler != nullptr);

		if (numVertices == 0) return;

		// prepare data for sampling a point on camera
		Ray2 ray;
		Spectrum contrib(1.0);

		// sample camera vertex!
		shared_ptr<CameraVertex> cameraVertex = make_shared<CameraVertex>();

		// sample camera position
		if (mDoStratification)
		{
			contrib = camera.samplePerPixelWe0ByPdf(cameraVertex.get(), nullptr, pixel, sampler->get1d());
		}
		else
		{
			contrib = camera.sampleWe0ByPdf(cameraVertex.get(), raster, nullptr, sampler->get1d());
		}
		ray.origin = cameraVertex->mPosition2;
		if (epCallbackFunc) epCallbackFunc(cameraVertex, contrib, 0.0);

		// sample camera direction
		contrib = camera.sampleWe1CosByPdf(&ray.direction, nullptr, *cameraVertex, sampler->get1d());
		assert(Math::IsNormalized(ray.direction));

		RandomWalk(sampler, *mScene, cameraVertex, ray, contrib, 0.0, numVertices - 1, callbackFunc);
	}
	
	shared_ptr<const Scene> mScene;
	bool mDoStratification;
};
