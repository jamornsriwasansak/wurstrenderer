#pragma once

#include "common/camera.h"
#include "common/scene.h"
#include "common/path.h"
#include "common/subpathsampler.h"
#include "common/sampler.h"

#include "common/util/file.h"

struct GeneralSubpathSampler : public SubpathSampler
{
	void randomWalk(Sampler * sampler,
					shared_ptr<Vertex> prevVertex,
					const Medium * medium,
					Ray3 ray,
					Spectrum prevContrib,
					double prevPdfW,
					const int numVertices,
					const Request requestFlags,
					const std::function<bool(shared_ptr<Vertex> & currentVertex,
											 shared_ptr<Vertex> & prevVertex,
											 const Vec3 & dirToCurrentVertex,
											 const double pdf)> & callbackFunc) const
	{
		assert(sampler != nullptr);

		if (numVertices == 0) return;

		int countVertices = 0;

		const bool needPdf = (requestFlags & Request::PdfA) != 0;

		while (true)
		{
			// trace ray to create next vertex
			shared_ptr<Vertex> vertex = mScene->intersect(&ray);

			// trace ray in medium
			if (mDoVolumetric)
			{
				Vec2 rndVec2 = sampler->get2d();
				if (medium)
				{
					// sample a point in volume
					Spectrum trContrib(1.0);
					shared_ptr<Vertex> mediumVertex = medium->sampleMediumVertex(&trContrib, ray, rndVec2);
					prevContrib *= trContrib;

					// assign vertex with medium vertex if distance sampling doesn't hit the surface
					if (mediumVertex) vertex = mediumVertex;
				}
			}

			// stop the path sampling process if a vertex can't be sampled.
			if (vertex == nullptr)
			{
				// callback
				if (callbackFunc) callbackFunc(vertex, prevVertex, -ray.direction, 0.0);
				break;
			}

			// assign vertex information
			vertex->mSubpathSampler = this;
			vertex->mPathContrib = prevContrib * prevVertex->mPathContrib;

			// compute pdf and convert to area domain if needed
			double pdf = prevPdfW;
			if (needPdf)
            {
			    pdf *= Vertex::ConvertToAreaFromSolidAngle(*prevVertex, *vertex);
                assert(pdf >= 0.0);
            }

			// callback
			if (callbackFunc) if(!callbackFunc(vertex, prevVertex, -ray.direction, pdf)) break;

			// num vertices exceeds, break!
			if (++countVertices == numVertices) break;

			if (vertex->hasType(VertexType::Surface)) // Hit surface
			{
				// init data for sampling next direction
				shared_ptr<SurfaceVertex> sv = Vertex::SafeCast<SurfaceVertex>(vertex);
				Vec3 inLocal = sv->mCoordFrame.toLocal(-ray.direction);
				Vec3 outLocal;

				// sample next direction 
				if (needPdf)
				{
					double pdfW;
					Spectrum bsdfCosByPdf = sv->mBsdf->sampleBsdfCosByPdf(&outLocal, &pdfW, *sv, inLocal, sampler->get2d());
					if (Math::IsZero(bsdfCosByPdf)) { break; }
					assert(Math::IsApprox(pdfW, sv->mBsdf->evalPdfW(*sv, inLocal, outLocal)));
					assert(Math::IsNormalized(outLocal));
					prevContrib = bsdfCosByPdf;
					prevPdfW = pdfW;
				}
				else
				{
					Spectrum bsdfCosByPdf = sv->mBsdf->sampleBsdfCosByPdf(&outLocal, nullptr, *sv, inLocal, sampler->get2d());
					if (Math::IsZero(bsdfCosByPdf)) { break; }
					assert(Math::IsNormalized(outLocal));
					prevContrib = bsdfCosByPdf;
				}

				// handle volume
				if (mDoVolumetric)
				{
					// define in and out medium
					sv->mMediumFrontface = (sv->mGeometry->mMediumFrontface == nullptr) ? medium : sv->mGeometry->mMediumFrontface.get();
					sv->mMediumBackface = (sv->mGeometry->mMediumBackface == nullptr) ? medium : sv->mGeometry->mMediumBackface.get();

					// assign medium just in case since the surface front and back medium interface could be set to nullptr
					if (Local3::IsFrontface(inLocal))
					{
						assert (sv->mGeometry->mMediumFrontface == nullptr || sv->mGeometry->mMediumFrontface.get() == medium);
						sv->mMediumFrontface = medium;
					}
					else
					{
						assert (sv->mGeometry->mMediumBackface == nullptr || sv->mGeometry->mMediumBackface.get() == medium);
						sv->mMediumBackface = medium;
					}
				}

				// update ray
				ray = Ray3(sv->mPosition, sv->mCoordFrame.toWorld(outLocal));
				prevVertex = vertex;
			}
			else if (vertex->hasType(VertexType::Medium))
			{
				assert(mDoVolumetric);
				shared_ptr<MediumVertex> mv = Vertex::SafeCast<MediumVertex>(vertex);
				Vec3 inLocal = mv->mCoordFrame.toLocal(-ray.direction);
				Vec3 outLocal;

				// sample next direction
				double pdfW;
				Spectrum phaseByPdf = mv->mPhaseFunction->samplePhaseByPdf(&outLocal, &pdfW, *mv, inLocal, sampler->get2d());
				if (Math::IsZero(phaseByPdf)) { break; }
				//assert(Math::IsApprox(pdfW, sv->mBsdf->evalPdfW(*sv, inLocal, outLocal)));
				assert(Math::IsNormalized(outLocal));
				prevContrib = phaseByPdf;
				prevPdfW = pdfW;

				// assign medium in and out
				mv->mMediumFrontface = medium;
				mv->mMediumBackface = medium;

				// update ray
				ray = Ray3(mv->mPosition, mv->mCoordFrame.toWorld(outLocal));
				prevVertex = vertex;
			}
			else // Hit light
			{
				break;
			}
		}
	}
	
	GeneralSubpathSampler(const shared_ptr<const Scene> & scene, const bool doCameraStratification = true, const bool doVolumetric = false):
		mScene(scene),
		mDoCameraStrafication(doCameraStratification),
		mDoVolumetric(doVolumetric)
	{
		// build cdftable weight
		std::vector<double> weights(mScene->mLightGeometries.size());
		for (Uint i = 0; i < mScene->mLightGeometries.size(); i++)
		{
			weights[i] = mScene->mLightGeometries[i]->mAreaLight->cdfWeight();
			assert(std::isfinite(weights[i]));
		}
		if (mScene->mEnvmapLight) weights.push_back(mScene->mEnvmapLight->cdfWeight());
		mLightCdfTable = CdfTable(weights);

		// map light pointer to cdftable index
		for (Uint iLight = 0; iLight < mScene->mLightGeometries.size(); iLight++)
		{
			mLightGeometryToCdfTableIndex[mScene->mLightGeometries[iLight].get()] = iLight;
		}
		if (mScene->mEnvmapLight) mLightGeometryToCdfTableIndex[nullptr] = mScene->mLightGeometries.size();
	}

	virtual void createLightSubpath(Sampler * sampler,
									const int numVertices,
									const Request requestFlags,
									const std::function<bool(shared_ptr<LightVertex> & lightVertex,
															 const double pdfA)> & epCallbackFunc,
									const std::function<bool(shared_ptr<Vertex> & currentVertex,
															 shared_ptr<Vertex> & prevVertex,
															 const Vec3 & dirToParent,
															 const double pdf)> & callbackFunc = nullptr) const
	{
		assert(sampler != nullptr);

		if (numVertices == 0) return;

		const bool needPdf = (requestFlags & Request::PdfA) != 0;

		// choose a light source
		double pdfLight;
		Vec2 sample = sampler->get2d();
		const Geometry * lightGeometry = sampleLight(&pdfLight, &sample[0], sample[0]);

		if (lightGeometry == nullptr) // envmap case
		{
			const Light * envmapLight = mScene->mEnvmapLight.get();
			shared_ptr<EnvmapVertex> envmapVertex = make_shared<EnvmapVertex>();
			envmapVertex->mGeometry = nullptr;

			// sample Le0 (outgoing direction of the enviroment map)
			Ray3 ray;
			double pdfDirection;

			const Spectrum contrib0 = envmapLight->sampleLe0ByPdf(envmapVertex.get(), &pdfDirection, sample) / pdfLight;
			envmapVertex->mSubpathSampler = this;
			envmapVertex->mPathContrib = contrib0;
			const double pdfLe0 = pdfLight * pdfDirection;

			// TODO:: handle volumetric case
			shared_ptr<LightVertex> lightVertex = envmapVertex;
			if (epCallbackFunc) if (!epCallbackFunc(lightVertex, pdfLe0)) return;

			// sample a point on the disk
			Vec3 diskPosition;
			double pdfDiskA;
			const Spectrum contrib1 = envmapLight->sampleLe1CosByPdf(&diskPosition, &pdfDiskA, *envmapVertex, sampler->get2d());
			ray.origin = diskPosition;
			ray.direction = -envmapVertex->mPosition;

			// continue the rest of the path
			randomWalk(sampler, lightVertex, nullptr, ray, contrib1, pdfDiskA, numVertices - 1, requestFlags, callbackFunc);
		}
		else // just any other light with geometry
		{
			const Light * lightsource = lightGeometry->mAreaLight.get();
			shared_ptr<LightVertex> lightVertex = make_shared<LightVertex>();
			lightVertex->mGeometry = lightGeometry;

			// sample Le0 (position on the lightsource)
			Ray3 ray;
			double pdfPosition;

			const Spectrum contrib0 = lightsource->sampleLe0ByPdf(lightVertex.get(), &pdfPosition, sample) / pdfLight;
			ray.origin = lightVertex->mPosition;
			lightVertex->mSubpathSampler = this;
			lightVertex->mPathContrib = contrib0;
			if (mDoVolumetric) lightVertex->mMediumBackface = lightsource->mMedium.get();
			const double pdfLe0 = pdfLight * pdfPosition;

			// callback
			if (epCallbackFunc) if (!epCallbackFunc(lightVertex, pdfLe0)) return;

			if (numVertices == 1) return;

			// sample Le1 (direction on from the lightsource)
			double pdfLe1 = 0.0;
			Vec3 outLocal;
			const Spectrum contrib1 = lightsource->sampleLe1CosByPdf(&outLocal, &pdfLe1, *lightVertex, sampler->get2d());
			ray.direction = lightVertex->mCoordFrame.toWorld(outLocal);
			assert(Math::IsNormalized(ray.direction));
			assert(Math::IsApprox(pdfLe1, lightsource->evalPdfLe1(*lightVertex, outLocal)));

			// continue the rest of the path
			randomWalk(sampler, lightVertex, lightsource->mMedium.get(), ray, contrib1, pdfLe1, numVertices - 1, requestFlags, callbackFunc);
		}
	}

	virtual void createCameraSubpath(Sampler * sampler,
									 Vec2 * raster,
									 const Camera & camera,
									 const Ivec2 & pixel,
									 const int numVertices,
									 const Request requestFlags,
									 const std::function<bool(shared_ptr<CameraVertex> & cameraVertex,
															  const double pdfA)> & epCallbackFunc,
									 const std::function<bool(shared_ptr<Vertex> & currentVertex,
															  shared_ptr<Vertex> & prevVertex,
															  const Vec3 & dirToParent,
															  const double pdfA)> & callbackFunc = nullptr) const
	{
		assert(sampler != nullptr);

		if (numVertices == 0) return;

		const bool needPdf = (requestFlags & Request::PdfA) != 0;

		// prepare data for sampling a point on camera
		Ray3 ray;
		double pdfAPos = 0.0, pdfWDir = 0.0;
		double * pdfAPtr = (needPdf) ? &pdfAPos : nullptr, * pdfWPtr = (needPdf) ? &pdfWDir : nullptr;

		// sample camera vertex!
		shared_ptr<CameraVertex> cameraVertex = make_shared<CameraVertex>();
		
		// sample camera position
		Spectrum contrib0 = camera.sampleWe0ByPdf(cameraVertex.get(), pdfAPtr, sampler->get2d());
		ray.origin = cameraVertex->mPosition;
		cameraVertex->mSubpathSampler = this;
		cameraVertex->mPathContrib = contrib0;
		if (mDoVolumetric) cameraVertex->mMediumBackface = camera.mMedium.get();

		// callback
		if (epCallbackFunc) if (!epCallbackFunc(cameraVertex, pdfAPos)) return;

		// sample camera direction
		Spectrum contrib1;
		if (mDoCameraStrafication)
		{
			contrib1 = camera.samplePerPixelWe1CosByPdf(&ray.direction, raster, pdfWPtr, *cameraVertex, pixel, sampler->get2d());
			if (pdfWPtr) assert(Math::IsApprox(*pdfWPtr, camera.evalPerPixelWe1PdfW(*cameraVertex, ray.direction)));
		}
		else
		{
			contrib1 = camera.sampleWe1CosByPdf(&ray.direction, raster, pdfWPtr, *cameraVertex, sampler->get2d());
			if (pdfWPtr) assert(Math::IsApprox(*pdfWPtr, camera.evalWe1PdfW(*cameraVertex, ray.direction)));
		}
		assert(Math::IsNormalized(ray.direction));

		// continue the rest of the path
		randomWalk(sampler, cameraVertex, camera.mMedium.get(), ray, contrib1, pdfWDir, numVertices - 1, requestFlags, callbackFunc);
	}

	double evalCameraPdfW(const CameraVertex & c, const Vec3 & w) const
	{
		return mDoCameraStrafication ? c.mCamera->evalPerPixelWe1PdfW(c, w) : c.mCamera->evalWe1PdfW(c, w);
	}

	double evalLightPdfA(const LightVertex & lv) const
	{
		// find index of light
		auto findResult = mLightGeometryToCdfTableIndex.find(lv.mGeometry);
		assert(findResult != mLightGeometryToCdfTableIndex.end());
		const Uint index = findResult->second;
		return mLightCdfTable.mPdfs[index] * lv.mLight->evalPdfLe0(lv);
	}

	const Geometry * sampleLight(double * pdfLight, double * remappedSample, const double sample) const
	{
		int lightIndex = mLightCdfTable.sample(pdfLight, remappedSample, sample);
		if (lightIndex < mScene->mLightGeometries.size()) return mScene->mLightGeometries[lightIndex].get();
		assert(lightIndex == mScene->mLightGeometries.size());
		return nullptr;
	}

	shared_ptr<const Scene> mScene;
	bool mDoVolumetric;
	bool mDoCameraStrafication;
	std::map<const Geometry *, Uint> mLightGeometryToCdfTableIndex;
	CdfTable mLightCdfTable;
};
