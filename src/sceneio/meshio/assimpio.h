#pragma once

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "common/wurst.h"
#include "common/bsdf.h"
#include "common/geometry.h"

#include "bsdf/lambertdiffuse.h"
#include "bsdf/phong.h"
#include "bsdf/mixed.h"
#include "bsdf/simplemetal.h"
#include "bsdf/twosided.h"
#include "camera/thinlens.h"
#include "geometry/trianglemesh.h"
#include "light/diffusearealight.h"
#include "light/smarealight.h"
#include "sampler/randomsampler.h"
#include "texture/constanttexture.h"
#include "texture/imagetexture.h"

#include "renderer/pathnextevent.h"
#include "renderer/bidirpath.h"
#include "renderer/matrixbidirpath.h"

struct AssimpInOut
{
	struct Result
	{
		std::vector<shared_ptr<TriangleMeshGeometry>> mGeometries;
		std::vector<shared_ptr<TriangleMeshGeometry>> mLightGeometries;
	};

	struct ObjTriMesh
	{
		shared_ptr<TriangleMeshGeometry> mTriMesh;
		int mMatId;
	};

	struct ObjBsdfTexture
	{
		shared_ptr<Texture<Spectrum>> mDiffuseReflectanceTexture;
		shared_ptr<Texture<Spectrum>> mSpecularReflectanceTexture;
		shared_ptr<Texture<double>> mExponentTexture;
		shared_ptr<Texture<Spectrum>> mEmissiveTexture;
		shared_ptr<Texture<double>> mAlphaTexture;
	};

	// TODO:: move obj related functions outside the assimpio
	static void SaveObj(const std::string & filename, const std::vector<Vec3> & points)
	{
		std::ofstream fileWriter;
		fileWriter.open(filename, std::ios::binary);
		for (const Vec3 &point : points)
			fileWriter << "v " << point[0] << " " << point[1] << " " << point[2] << std::endl;
		fileWriter.close();
	}

	static Result LoadFromObj(const filesystem::path & objFilepath,
							  const bool doConvertSrgbTheMitsubaWay = false,
							  const shared_ptr<const Bsdf> overrideBsdf = nullptr,
							  const shared_ptr<Light> overrideArealight = nullptr,
							  const bool doAverageFaceNormals = false,
							  const int dim = 3,
							  const double planeZ = 0.0)
	{
		assert(dim == 3 || dim == 2);

		Result result;

		std::vector<ObjTriMesh> otms;
		std::vector<ObjBsdfTexture> obts;
		PopulateObjVectors(objFilepath, &otms, &obts, doConvertSrgbTheMitsubaWay, doAverageFaceNormals);

		std::vector<shared_ptr<const Bsdf>> bsdfs;
		std::vector<shared_ptr<Light>> areaLights;
		for (const ObjBsdfTexture & obt : obts)
		{
			// arealight
			if (!Math::IsZero(obt.mEmissiveTexture->average()))
			{
				if (obt.mExponentTexture->average() == 0.0)
				{
					areaLights.push_back(make_shared<DiffuseAreaLight>(obt.mEmissiveTexture));
				}
				else
				{
					areaLights.push_back(make_shared<SimpleMetalAreaLight>(obt.mEmissiveTexture, SimpleMetalBsdf::AngleFromPhongExponent(obt.mExponentTexture->average())));
				}
			}
			else
			{
				areaLights.push_back(nullptr);
			}

			shared_ptr<const Bsdf> diffuseBsdf = make_shared<const LambertDiffuseBsdf>(obt.mDiffuseReflectanceTexture);
			shared_ptr<const Bsdf> glossyBsdf = nullptr;

			if (dim == 3) glossyBsdf = make_shared<const PhongBsdf>(obt.mSpecularReflectanceTexture, obt.mExponentTexture);
			else if (dim == 2) glossyBsdf = make_shared<const SimpleMetalBsdf>(obt.mSpecularReflectanceTexture, obt.mExponentTexture);

			// material
			if (Math::IsZero(obt.mDiffuseReflectanceTexture->average())) // has specular component only
			{
				bsdfs.push_back(make_shared<const TwoSidedBsdf>(glossyBsdf));
			}
			else if (Math::IsZero(obt.mSpecularReflectanceTexture->average())) // has diffuse component only
			{
				bsdfs.push_back(make_shared<const TwoSidedBsdf>(diffuseBsdf));
			}
			else // has both component
			{
				std::vector<shared_ptr<const Bsdf>> twoBsdfs;
				twoBsdfs.emplace_back(diffuseBsdf);
				twoBsdfs.emplace_back(glossyBsdf);
				shared_ptr<const MixedBsdf> mixedBsdf = make_shared<const MixedBsdf>(twoBsdfs);
				bsdfs.push_back(make_shared<const TwoSidedBsdf>(mixedBsdf));
			}
		}

		for (ObjTriMesh & otm : otms)
		{
			shared_ptr<Light> areaLight = (overrideArealight == nullptr) ? areaLights[otm.mMatId] : overrideArealight ;
			shared_ptr<const Bsdf> bsdf = (overrideBsdf == nullptr) ? bsdfs[otm.mMatId] : overrideBsdf;

			if (areaLight != nullptr)
			{
				otm.mTriMesh->mAreaLight = areaLight;
				if (dim == 2)
                {
                    otm.mTriMesh->precomputeForSampling2(planeZ);
                }
				else if (dim == 3)
				{
					otm.mTriMesh->precomputeForSampling();
				}
				result.mLightGeometries.push_back(otm.mTriMesh);
			}
			else
            {
                if (bsdf != nullptr)
                {
                    otm.mTriMesh->mBsdf = bsdf;
                    otm.mTriMesh->mAlpha = obts[otm.mMatId].mAlphaTexture;
                }
                result.mGeometries.push_back(otm.mTriMesh);
            }
		}

		return result;
	}

	static void PopulateObjVectors(const filesystem::path & meshFilepath,
								   std::vector<ObjTriMesh> * triMesh,
								   std::vector<ObjBsdfTexture> * texture,
								   const bool doConvertSrgbTheMitsubaWay = false,
								   const bool doAverageFaceNormals = false)
	{
		unsigned int aiProcesses = aiProcessPreset_TargetRealtime_Fast;

		std::vector<shared_ptr<const Texture<double>>> alphaTextures;

		// load aiScene
		Assimp::Importer importer;
		const aiScene * scene = importer.ReadFile(meshFilepath.string().c_str(), aiProcesses);
		if (scene == nullptr) throw std::runtime_error("cannot open file");

		// translate textures
		const int numMaterials = scene->mNumMaterials;
		for (int iMat = 0; iMat < numMaterials; iMat++)
		{
			// note: assimp always generate material index = 0 as DefaultMaterial automatically
			const auto aiMat = scene->mMaterials[iMat];
			ObjBsdfTexture obt;

			AimatResult diffTexture = AimatToTexture(meshFilepath, *aiMat, aiTextureType_DIFFUSE, AI_MATKEY_COLOR_DIFFUSE, doConvertSrgbTheMitsubaWay);
			obt.mDiffuseReflectanceTexture = diffTexture.mSpectrum;
			obt.mSpecularReflectanceTexture = AimatToTexture(meshFilepath, *aiMat, aiTextureType_SPECULAR, AI_MATKEY_COLOR_SPECULAR, doConvertSrgbTheMitsubaWay).mSpectrum;
			obt.mExponentTexture = AimatToTexture(meshFilepath, *aiMat, aiTextureType_SHININESS, AI_MATKEY_SHININESS, false).mAlpha;
			obt.mEmissiveTexture = AimatToTexture(meshFilepath, *aiMat, aiTextureType_EMISSIVE, AI_MATKEY_COLOR_EMISSIVE, false).mSpectrum;
			obt.mAlphaTexture = diffTexture.mAlpha;
			texture->emplace_back(obt);
		}

		// translate trimesh
		for (int iMesh = 0; iMesh < int(scene->mNumMeshes); iMesh++)
		{
			ObjTriMesh otm;
			otm.mTriMesh = make_shared<TriangleMeshGeometry>();
			const auto aiSceneMesh = scene->mMeshes[iMesh];
			const int numVertices = aiSceneMesh->mNumVertices;
			const bool hasTextureCoords = aiSceneMesh->HasTextureCoords(0);
			const bool hasShadingNormal = aiSceneMesh->HasNormals();
			otm.mTriMesh->mPositions.reserve(numVertices);
			otm.mTriMesh->mNormals.reserve(numVertices);
			otm.mTriMesh->mTextureCoords.reserve(numVertices);

			for (int iVert = 0; iVert < numVertices; iVert++)
			{
				// position
				{
					const float x = aiSceneMesh->mVertices[iVert][0];
					const float y = aiSceneMesh->mVertices[iVert][1];
					const float z = aiSceneMesh->mVertices[iVert][2];
					otm.mTriMesh->mPositions.emplace_back(x, y, z);
				}

				// normal
				if (hasShadingNormal)
				{
					const float x = aiSceneMesh->mNormals[iVert][0];
					const float y = aiSceneMesh->mNormals[iVert][1];
					const float z = aiSceneMesh->mNormals[iVert][2];
					if (Math::IsApprox(Math::Length(Vec3(x, y, z)), 0.0))
					{
						otm.mTriMesh->mNormals.emplace_back(Vec3(0.0, 0.0, 0.0));
					}
					else
					{
						otm.mTriMesh->mNormals.emplace_back(Math::Normalize(Vec3(x, y, z)));
					}
				}

				// tex coords
				if (hasTextureCoords)
				{
					const float u = aiSceneMesh->mTextureCoords[0][iVert][0];
					const float v = aiSceneMesh->mTextureCoords[0][iVert][1];
					otm.mTriMesh->mTextureCoords.emplace_back(u, v);
				}
				else
				{
					otm.mTriMesh->mTextureCoords.emplace_back(0, 0);
				}
			}

			const int numTriangles = aiSceneMesh->mNumFaces;
			otm.mTriMesh->mTriangles.reserve(numTriangles);
			for (int iIdx = 0; iIdx < numTriangles; iIdx++)
			{
				unsigned int t0 = aiSceneMesh->mFaces[iIdx].mIndices[0];
				unsigned int t1 = aiSceneMesh->mFaces[iIdx].mIndices[1];
				unsigned int t2 = aiSceneMesh->mFaces[iIdx].mIndices[2];
				otm.mTriMesh->mTriangles.emplace_back(t0, t1, t2);
			}

			otm.mMatId = aiSceneMesh->mMaterialIndex;

			// override all loaded normals if smooth normals are required
			if (doAverageFaceNormals)
			{
				otm.mTriMesh->mNormals = std::vector<Vec3>(numVertices, Vec3(0.0));
				for (int iTriangle = 0; iTriangle < numTriangles; iTriangle++)
				{
					Ivec3 indices = otm.mTriMesh->mTriangles[iTriangle];

					// fetch 3 positions of a triangle
					Vec3 positions[3];
					for (int i = 0; i < 3; i++) positions[i] = otm.mTriMesh->mPositions[indices[i]];

					// compute triangle normal
					Vec3 normal = Math::Normalize(Math::Cross(positions[1] - positions[0], positions[2] - positions[0]));

					// accumulate triangle normal to the vertices
					for (int i = 0; i < 3; i++) otm.mTriMesh->mNormals[indices[i]] += normal;
				}

				for (int iNormal = 0; iNormal < numVertices; iNormal++)
				{
					Vec3 & normal = otm.mTriMesh->mNormals[iNormal];
					normal = Math::IsApprox(Math::Length2(normal), 0.0) ? Vec3(0.0, 0.0, 0.0) : Math::Normalize(normal);
				}
			}

			triMesh->emplace_back(otm);
		}
	}

	struct AimatResult
	{
		shared_ptr<Texture<Spectrum>> mSpectrum = nullptr;
		shared_ptr<Texture<double>> mAlpha = nullptr;
	};

	static AimatResult AimatToTexture(const filesystem::path & meshFilepath,
									  const aiMaterial & aiMat,
									  const aiTextureType & textureKey,
									  const char * colorKey,
									  unsigned int type,
									  unsigned int index,
									  const bool doConvertSrgbTheMitsubaWay)
	{
		aiString aiTexturePath;
		aiColor3D aiReflectance;

		AimatResult result;

		if (aiMat.Get(AI_MATKEY_TEXTURE(textureKey, 0), aiTexturePath) == aiReturn_SUCCESS)
		{
			filesystem::path textureFilepath = meshFilepath.parent_path() / filesystem::path(aiTexturePath.C_Str());
			int numChannels = 0;
			result.mSpectrum = make_shared<ImageTexture<Spectrum>>(FimageIo::LoadPtr(&numChannels, textureFilepath));
			result.mAlpha = (numChannels == 4) ? make_shared<ImageTexture<double>>(FimageIo::LoadPtr(nullptr, textureFilepath, 3)) : nullptr;
			return result;
		}
		else if (aiMat.Get(colorKey, (unsigned int)type, (unsigned int)index, aiReflectance) == aiReturn_SUCCESS)
		{
			if (textureKey == aiTextureType::aiTextureType_SHININESS)
			{
				// fix shininess bug multiplied by 4 introduced by assimp (they said "to match what most renderers do")
				double color = double(aiReflectance.r);
				result.mAlpha = make_shared<ConstantTexture<double>>(color / 4.0);
			}
			else
			{
				Spectrum color;
				if (doConvertSrgbTheMitsubaWay)
				{
					color = SpectrumConvertUtil::SpectrumFromSrgb(static_cast<double>(aiReflectance.r),
																  static_cast<double>(aiReflectance.g),
																  static_cast<double>(aiReflectance.b));
				}
				else
				{
					color = SpectrumConvertUtil::SpectrumFromRgb(static_cast<double>(aiReflectance.r),
																 static_cast<double>(aiReflectance.g),
																 static_cast<double>(aiReflectance.b));
				}
				result.mSpectrum = make_shared<ConstantTexture<Spectrum>>(color);
			}
			return result;
		}

		throw std::runtime_error("loading texture error");
	}
};
