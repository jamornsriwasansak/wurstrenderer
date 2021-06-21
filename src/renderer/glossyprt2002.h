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

struct GlossyPrt2002 : public SingleCameraRenderer
{
	static shared_ptr<GlossyPrt2002> FromSimple2Json(shared_ptr<const Camera> & camera,
												 shared_ptr<const Scene> & scene,
												 shared_ptr<Sampler> & sampler,
												 const json & json,
												 const filesystem::path & rootJsonPath)
	{
		const int numAaSamples = Util::Json::GetValue<int>(json, "num_aa_samples", 16);
		const int numSamplesPerPoint = Util::Json::GetValue<int>(json, "num_samples_per_point", 10000);
		const int numEnvmapCoeffs = Util::Json::GetValue<int>(json, "num_envmap_bands", 3);
		const int numLocalCoeffs = Util::Json::GetValue<int>(json, "num_local_bands", 3);
		const int numPathVertices = Util::Json::GetValue<int>(json, "num_path_vertices", 3);
		return make_shared<GlossyPrt2002>(camera, scene, sampler, numSamplesPerPoint, numEnvmapCoeffs, numLocalCoeffs, numAaSamples, numPathVertices);
	}

	GlossyPrt2002(shared_ptr<const Camera> & camera,
				  shared_ptr<const Scene> & scene,
				  shared_ptr<Sampler> & sampler,
				  const int numSamplesPerPoint,
				  const int numEnvmapShBand,
				  const int numLocalShBand,
				  const int numAaSamples,
				  const int numPathVertices) :
		SingleCameraRenderer(camera, scene, sampler),
		mVerticesVisitor(scene),
		mNumEnvmapCoeffs(Math::Square(numEnvmapShBand)),
		mNumLocalCoeffs(Math::Square(numLocalShBand)),
		mNumMatrixElements(Math::Square(numEnvmapShBand * numLocalShBand)),
		mNumSamplesPerPoint(numSamplesPerPoint),
		mNumAaSamples(numAaSamples),
		mNumPathVertices(numPathVertices)
	{
	}

	int64_t globalIndex(const int iMesh, const int iPoint) const
	{
		const int globaIVertex = mVerticesVisitor.mCumulativeNumPoints[iMesh] + iPoint;
		return int64_t(globaIVertex) * int64_t(mNumMatrixElements);
	}

	void precomputeTransferMatrices() const
	{
		GeneralSubpathSampler pathSampler(mScene, false, false);
		double invNumSamplesPerPoint = 1.0 / double(mNumSamplesPerPoint);

		int64_t numElements = int64_t(mVerticesVisitor.mNumTotalPoints) * int64_t(mNumMatrixElements);
		mTransferMatrices = std::vector<Spectrum>(numElements, Spectrum(0.0));

		shared_ptr<ProgressReport> progress = make_shared<ProgressReport>(mVerticesVisitor.mNumTotalPoints);
		Viewer::SetProgressReport(mCamera->mFilm, progress);

		mVerticesVisitor.forEachPoint<Parallel>([&](const int iMesh, shared_ptr<const TriangleMeshGeometry> & triMesh, const int iPoint)
		{
			// get sampler and surface vertex
			shared_ptr<Sampler> sampler = mSampler->clone(int(globalIndex(iMesh, iPoint)));
			shared_ptr<SurfaceVertex> sv = mVerticesVisitor.getSurfaceVertex(iMesh, iPoint);

			int64_t iMatrixBegin = globalIndex(iMesh, iPoint);

			std::vector<double> localShValues(mNumLocalCoeffs);
			std::vector<double> envmapShValues(mNumEnvmapCoeffs);

			for (int i = 0; i < mNumSamplesPerPoint; i++)
			{
				// sampling incoming radiance to the surface
				Vec3 sampled = Mapping::HemisphereFromSquare(sampler->get2d());
				Vec3 initialLocalDirection = sv->mCoordFrame.toWorld(sampled);
				double initialPdf = 0.5 * Math::InvPi;
				Ray out(sv->mPosition, initialLocalDirection);

				pathSampler.randomWalk(sampler.get(), sv, nullptr, out, Spectrum(1.0 / initialPdf), 0.0, mNumPathVertices - 2, SubpathSampler::Request::None,
					[&](shared_ptr<Vertex> & currentVertex, shared_ptr<Vertex> & prevVertex, const Vec3 &, const double pdf) -> bool
					{
						assert(currentVertex != nullptr);
						if (currentVertex->hasType(VertexType::Light))
						{
							// out going direction to envmap
							Vec3 envmapDirection = Math::Normalize(currentVertex->mPosition - prevVertex->mPosition);

							if (mNumPathVertices == 3) // this means direct light
							{
								assert(Math::IsApprox(Math::Dot(envmapDirection, initialLocalDirection), 1.0));
							}

							// evaluate all spherical harmonics terms for local
							for (int iCoeff = 0; iCoeff < mNumLocalCoeffs; iCoeff++)
							{
								localShValues[iCoeff] = Math::SphericalHarmonics::Eval(iCoeff, initialLocalDirection);
							}

							// evaluate all spherical harmonics terms for envmap
							for (int iCoeff = 0; iCoeff < mNumEnvmapCoeffs; iCoeff++)
							{
								envmapShValues[iCoeff] = Math::SphericalHarmonics::Eval(iCoeff, envmapDirection);
							}

							for (int iMatrixElement = 0; iMatrixElement < mNumMatrixElements; iMatrixElement++)
							{
								int row = iMatrixElement / mNumLocalCoeffs;
								int col = iMatrixElement % mNumLocalCoeffs;
								mTransferMatrices[iMatrixBegin + iMatrixElement] += localShValues[col] * envmapShValues[row] * currentVertex->mPathContrib * invNumSamplesPerPoint;
							}
						}
						return true;
					});
			}

			// update percentage on the viewer
			progress->increment(1);
			Viewer::Redraw(mCamera->mFilm, Ibound2(Ivec2(0, 0), Ivec2(0, 0)));
		});
		Util::Vector::WriteVectorBin("test.gprt2002.bin", mTransferMatrices);
	}

	void interpolateMatrix(std::vector<Spectrum> * matrix, std::vector<Spectrum> & matrices, const TriangleInfoVertex & tiv) const
	{
		int64_t begin[3];

		// find begin and end
		for (int i = 0; i < 3; i++)
		{
			begin[i] = globalIndex(tiv.mGeometryId, tiv.mIndices[i]);
		}

		for (int i = 0; i < mNumMatrixElements; i++)
		{
			Spectrum element[3];
			for (int j = 0; j < 3; j++)
			{
				element[j] = matrices[begin[j] + i];
			}
			(*matrix)[i] = Math::BarycentricInterpolate(element[0], element[1], element[2], tiv.mBarycentricCoords);
		}
	}

	void matrixVectorMultiply(std::vector<Spectrum> * result, const std::vector<Spectrum> & mat, const std::vector<Spectrum> vec) const
	{
		for (int r = 0; r < mNumEnvmapCoeffs; r++)
		{
			Spectrum x(0.0);
			for (int c = 0; c < mNumLocalCoeffs; c++)
			{
				x += mat[r * mNumLocalCoeffs + c] * vec[c];
			}
			(*result)[r] = x;
		}
	}

	std::vector<Spectrum> precomputeEnvmap() const
	{
		shared_ptr<const EquiRectRep<Spectrum>> equirectEnvmap = static_pointer_cast<const EquiRectRep<Spectrum>>(mScene->mEnvmapLight->mSphRep);
		FimageIo::Save(SphRepConvertUtil::EquiRectRepFrom(SphRepConvertUtil::SphHarRepFrom(*equirectEnvmap, mNumEnvmapCoeffs), Ivec2(1024, 512)).mImage, "checkcheck.pfm");
		return SphRepConvertUtil::SphHarRepFrom(*equirectEnvmap, mNumEnvmapCoeffs).mCoeffs;
	}

	void execute(const std::vector<Spectrum> & envmapCoeffs) const
	{
		mCamera->mFilm->requestBaseBuffer();

		// act as a rasterizer
		Parallel::Split2d(mCamera->mFilm->mResolution, [&](const Ibound2 & bound)
		{
			int tileSeed = bound.pMin[0] * mCamera->mFilm->mResolution[1] + bound.pMin[1];
			shared_ptr<Sampler> tileSampler = mSampler->clone(tileSeed);

			std::vector<Spectrum> interpolatedMatrix(mNumMatrixElements);

			std::vector<Spectrum> localCoeffs(mNumLocalCoeffs);

			for (const Ivec2 & pixel : bound)
			{
				Spectrum result(0.0);
				for (Uint iSample = 0; iSample < mNumAaSamples; iSample++)
				{
					// sample camera ray
					Ray3 r;
					CameraVertex cv;
					mCamera->sampleWe0ByPdf(&cv, nullptr, tileSampler->get2d());
					r.origin = cv.mPosition;
					mCamera->samplePerPixelWe1CosByPdf(&r.direction, nullptr, nullptr, cv, pixel, tileSampler->get2d());
					
					// intersect triangle
					TriangleInfoVertex vertex;

					Ray3 r1 = r;
					Ray3 r2 = r;
					if (mScene->mAccel->intersectTriangleInfo(&vertex, &r1))
					{
						// get matrix at the shading point
						interpolateMatrix(&interpolatedMatrix, mTransferMatrices, vertex);

						// compute transfer
						matrixVectorMultiply(&localCoeffs, interpolatedMatrix, envmapCoeffs);

						// debug
						shared_ptr<SurfaceVertex> v = Vertex::SafeCast<SurfaceVertex>(mScene->intersect(&r2));
						Vec3 in = v->mCoordFrame.toLocal(-r.direction);
						Vec3 out;
						Spectrum bsdfContrib = v->mBsdf->sampleBsdfCosByPdf(&out, nullptr, *v, in, tileSampler->get2d());
						out = v->mCoordFrame.toWorld(out);

						// reconstruct radiance
						for (int iCoeff = 0; iCoeff < mNumLocalCoeffs; iCoeff++)
						{
							const double shVal = Math::SphericalHarmonics::Eval(iCoeff, out);
							for (int iChannel = 0; iChannel < Spectrum::NumElements; iChannel++)
							{
								result[iChannel] += shVal * bsdfContrib[iChannel] * localCoeffs[iCoeff][iChannel];
							}
						}
					}
					else
					{
						// reconstruct envmap from sh values
						for (int iCoeff = 0; iCoeff < mNumEnvmapCoeffs; iCoeff++)
						{
							const double shVal = Math::SphericalHarmonics::Eval(iCoeff, -r.direction);
							for (int iChannel = 0; iChannel < Spectrum::NumElements; iChannel++)
							{
								result += shVal * envmapCoeffs[iCoeff][iChannel];
							}
						}
					}
				}
				mCamera->mFilm->addSample(pixel, Spectrum(result) / static_cast<double>(mNumAaSamples));
			}
			Viewer::Redraw(mCamera->mFilm, bound);
		});
	}

	void render() override
	{
		precomputeTransferMatrices();
		filesystem::path k = "test.gprt2002.bin";
		mTransferMatrices = Util::Vector::LoadVectorBin<Spectrum>(k);
		std::vector<Spectrum> envmapCoeffs = precomputeEnvmap();
		execute(envmapCoeffs);
	}

	int mNumEnvmapCoeffs;
	int mNumLocalCoeffs;
	int mNumMatrixElements;
	int mNumSamplesPerPoint;
	int mNumAaSamples;
	int mNumPathVertices;
	TriangleMeshVerticesVisitor mVerticesVisitor;

	mutable std::vector<Spectrum> mTransferMatrices;
};