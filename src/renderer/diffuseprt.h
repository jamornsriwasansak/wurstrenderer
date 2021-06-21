#pragma once

#include "common/wurst.h"

#include "common/renderer.h"
#include "common/scene.h"
#include "common/trimeshvertexiterator.h"
#include "common/util/vector.h"
#include "common/util/json.h"
#include "common/viewer.h"

#include "pathsampler/generalsubpathsampler.h"
#include "sphericalrepresentation/sphericalharmonics.h"
#include "sphericalrepresentation/equirectangular.h"
#include "sphericalrepresentation/convertutil.h"

struct DiffusePrt : public SingleCameraRenderer
{
	static shared_ptr<DiffusePrt> FromSimple2Json(shared_ptr<const Camera> & camera,
												  shared_ptr<const Scene> & scene,
												  shared_ptr<Sampler> & sampler,
												  const json & json,
												  const filesystem::path & rootJsonPath)
	{
		const int numPathVertices = Util::Json::GetValue<int>(json, "num_path_vertices", 10);
		const int numAaSamples = Util::Json::GetValue<int>(json, "num_aa_samples", 16);
		const int numSamplesPerPoint = Util::Json::GetValue<int>(json, "num_samples_per_point", 1000);
		const int numCoeffs = Util::Json::GetValue<int>(json, "num_coeffs", 8);
		return make_shared<DiffusePrt>(camera, scene, sampler, numPathVertices, numSamplesPerPoint, numCoeffs, numAaSamples);
	}

	DiffusePrt(shared_ptr<const Camera> & camera,
			   shared_ptr<const Scene> & scene,
			   shared_ptr<Sampler> & sampler,
			   const int numVertices,
			   const int numSamplesPerBakingPoint,
			   const int numCoeffs,
			   const int numAaSamples):
		SingleCameraRenderer(camera, scene, sampler),
		mNumPathVertices(numVertices),
		mNumSamplesPerPoint(numSamplesPerBakingPoint),
		mNumCoeffs(numCoeffs),
		mNumAaSamples(numAaSamples),
		mVerticesVisitor(scene)
	{
		// check if the scene file is contaminated with weird envmap
		if (scene->mLightGeometries.size() != 0)
		{
			throw std::runtime_error("prt should contain only one lightsource which is envmap");
		}

		if (scene->mEnvmapLight == nullptr)
		{
			throw std::runtime_error("prt should have an environment map");
		}
	}

	int globalCoeffIndex(const int iMesh, const int iPoint, const int iCoeff) const
	{
		const int globalIVertex = mVerticesVisitor.mCumulativeNumPoints[iMesh] + iPoint;
		return (globalIVertex * mNumCoeffs) + iCoeff;
	}

	void reverseGlobalIndex(int * iMesh, int * iPoint, int * iCoeff, const int globalIndex) const
	{
		int globalIVertex = globalIndex / mNumCoeffs;
		auto iter = std::upper_bound(mVerticesVisitor.mCumulativeNumPoints.begin(), mVerticesVisitor.mCumulativeNumPoints.end(), globalIVertex);
		--iter;
		*iMesh = int(std::distance(mVerticesVisitor.mCumulativeNumPoints.begin(), iter));
		*iPoint = globalIVertex - *iter;
		*iCoeff = globalIndex - globalIVertex * mNumCoeffs;
	}

	void precomputeForGeometries() const
	{
		GeneralSubpathSampler pathSampler(mScene, false, false);
		double invNumSamplesPerPoint = 1.0 / double(mNumSamplesPerPoint);

		mGeometriesCoeffs = std::vector<Spectrum>(mVerticesVisitor.mNumTotalPoints * mNumCoeffs, Spectrum(0.0));

		shared_ptr<ProgressReport> progress = make_shared<ProgressReport>(mVerticesVisitor.mNumTotalPoints);
		Viewer::SetProgressReport(mCamera->mFilm, progress);
		mVerticesVisitor.forEachPoint<Parallel>([&](int iMesh, shared_ptr<const TriangleMeshGeometry> & triMesh, int iPoint)
		{
			shared_ptr<Sampler> sampler = mSampler->clone(0);

			Vec3 position = triMesh->mPositions[iPoint];
			Vec3 normal = triMesh->mNormals[iPoint];
			Vec2 texCoord = triMesh->mTextureCoords[iPoint];

			// initialize surface vertex
			shared_ptr<SurfaceVertex> sv = make_shared<SurfaceVertex>();
			sv->mBsdf = triMesh->mBsdf.get();
			sv->mCoordFrame = CoordFrame3(normal);
			sv->mGeometry = triMesh.get();
			sv->mGeometryNormal = normal;
			sv->mPathContrib = Spectrum(1.0);
			sv->mPosition = position;
			sv->mSubpathSampler = &pathSampler;
			sv->mTextureCoord = texCoord;

			for (int i = 0;i < mNumSamplesPerPoint;i++)
			{ 
				// sample a direction (this is our local outgoing radiance).
				Vec3 sampledDirection = Mapping::HemisphereFromSquare(sampler->get2d());
				Vec3 outgoingDirection;
				Spectrum initialContrib = sv->mBsdf->sampleBsdfCosByPdf(&outgoingDirection, nullptr, *sv, sampledDirection, sampler->get2d());
				Ray3 out(position, sv->mCoordFrame.toWorld(outgoingDirection));

				// random walk
				pathSampler.randomWalk(sampler.get(), sv, nullptr, out, initialContrib, 0.0, mNumPathVertices - 2, SubpathSampler::Request::None,
					[&](shared_ptr<Vertex> currentVertex, shared_ptr<Vertex> prevVertex, const Vec3 &, const double pdf) -> bool
					{
						assert(currentVertex != nullptr);
						if (currentVertex->hasType(VertexType::Light))
						{
							// outgoing direction to envmap
							Vec3 outgoing = Math::Normalize(currentVertex->mPosition - prevVertex->mPosition);

							// loop over and update all coeffs
							for (int iCoeff = 0; iCoeff < mNumCoeffs; iCoeff++)
							{
								// eval spherical harmonics for outgoing direction
								Ivec2 lm = Math::SphericalHarmonics::LmFromShIndex(iCoeff);
								double shVal = Math::SphericalHarmonics::Eval(lm[0], lm[1], outgoing);
								mGeometriesCoeffs[globalCoeffIndex(iMesh, iPoint, iCoeff)] += shVal * currentVertex->mPathContrib * invNumSamplesPerPoint;
							}
						}
						return true;
					});
			}

			// update percentage on the viewer
			progress->increment(1);
			Viewer::Redraw(mCamera->mFilm, Ibound2(Ivec2(0, 0), Ivec2(0, 0)));
		});

		// save coefficients
		Util::Vector::WriteVectorBin<Spectrum>("test.prt.bin", mGeometriesCoeffs);
	}

	void precomputeForEnvmap() const
	{
		shared_ptr<const EquiRectRep<Spectrum>> equirectEnvmap = static_pointer_cast<const EquiRectRep<Spectrum>>(mScene->mEnvmapLight->mSphRep);
		mEnvmapCoeffs = SphRepConvertUtil::SphHarRepFrom(*equirectEnvmap, mNumCoeffs).mCoeffs;
	}

	Spectrum getInterpolatedCoeffs(const TriangleInfoVertex & tiv, const int iCoeff) const
	{
		Spectrum triCoeff[3];
		for (int i = 0; i < 3; i++)
		{
			int index = globalCoeffIndex(tiv.mGeometryId, tiv.mIndices[i], iCoeff);
			triCoeff[i] = mGeometriesCoeffs[index];
		}
		return Math::BarycentricInterpolate(triCoeff[0], triCoeff[1], triCoeff[2], tiv.mBarycentricCoords);
	}

	void execute(const int coeffIndex, const bool visualize) const
	{
		mCamera->mFilm->requestBaseBuffer();

		mGeometriesCoeffs = Util::Vector::LoadVectorBin<Spectrum>("test.prt.bin");

		// act as rasterizer
		Parallel::Split2d(mCamera->mFilm->mResolution, [&](const Ibound2 & bound)
		{
			int tileSeed = bound.pMin[0] * mCamera->mFilm->mResolution[1] + bound.pMin[1];
			RandomSampler tileSampler(tileSeed);

			for (const Ivec2 & pixel : bound)
			{
				Spectrum result(0.0);
				for (Uint iSample = 0; iSample < mNumAaSamples; iSample++)
				{
					// sample camera ray
					Ray3 r;
					CameraVertex cv;
					mCamera->sampleWe0ByPdf(&cv, nullptr, tileSampler.get2d());
					r.origin = cv.mPosition;
					mCamera->samplePerPixelWe1CosByPdf(&r.direction, nullptr, nullptr, cv, pixel, tileSampler.get2d());
					
					// intersect triangle
					TriangleInfoVertex vertex;
					if (mScene->mAccel->intersectTriangleInfo(&vertex, &r))
					{
						if (visualize)
						{
							// reconstruct
							result += getInterpolatedCoeffs(vertex, coeffIndex);
						}
						else
						{
							for (int iCoeff = 0; iCoeff < mNumCoeffs; iCoeff++)
							{
								// Dot product is here
								result += getInterpolatedCoeffs(vertex, iCoeff) * mEnvmapCoeffs[iCoeff];
							}
						}
					}
					else
					{
						// TODO:: render envmap
						result += 0;
					}
				}
				if (visualize)
				{
					result = (result[0] < 0.0) ? Spectrum(-result[0], 0.0, 0.0) : Spectrum(0.0, result[0], 0.0);
				}
				mCamera->mFilm->addSample(pixel, Spectrum(result) / static_cast<double>(mNumAaSamples));
			}
			Viewer::Redraw(mCamera->mFilm, bound);
		});
	}

	void render() override
	{
		precomputeForGeometries();
		precomputeForEnvmap();
		for (int i = 0; i < mNumCoeffs; i++)
		{
			execute(i, true);
			FimageIo::Save(mCamera->mFilm->getImage(), "vis" + std::to_string(i) + ".pfm");
			mCamera->mFilm->clear();
		}
		execute(0, false);
	}

	// initial parameters
	filesystem::path mShCoeffCacheFilepath;
	int mNumPathVertices;
	int mNumSamplesPerPoint;
	int mNumCoeffs;
	int mNumAaSamples;
	TriangleMeshVerticesVisitor mVerticesVisitor;

	mutable std::vector<Spectrum> mGeometriesCoeffs;
	mutable std::vector<Spectrum> mEnvmapCoeffs;
};
