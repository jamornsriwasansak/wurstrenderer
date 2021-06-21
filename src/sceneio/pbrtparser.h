#pragma once

#include "common/wurst.h"

#include "common/scene.h"
#include "common/camera.h"
#include "common/viewer.h"

#include "accel/embreeaccel.h"
#include "bsdf/lambertdiffuse.h"
#include "bsdf/orennayar.h"
#include "bsdf/phong.h"
#include "bsdf/twosided.h"
#include "camera/thinlens.h"
#include "geometry/cube.h"
#include "geometry/plane.h"
#include "geometry/trianglemesh.h"
#include "light/diffusearealight.h"
#include "light/smarealight.h"
#include "renderer/pathnextevent.h"
#include "renderer/lightsplatting.h"
#include "renderer/photonmap.h"
#include "renderer/primitivepath.h"
#include "renderer/pssmltprimitivepath.h"
#include "renderer/visualizegeomnormal.h"
#include "sampler/randomsampler.h"
#include "sphericalrepresentation/equirectangular.h"
#include "sceneio/meshio/assimpio.h"
#include "sceneio/meshio/wo3io.h"
#include "texture/constanttexture.h"
#include "texture/checkerboard.h"
#include "texture/imagetexture.h"

#include "pbrtParser/Scene.h"

// Ingo Wald Pbrt Parser
using std::cout;
using std::endl;
// the output file we're writing.
std::ofstream out;
// FILE *out = NULL;

size_t numWritten = 0;
size_t numVerticesWritten = 0;

struct PbrtParser
{
	static shared_ptr<Texture<Spectrum>> ParseTexture(const pbrt::vec3f & pbrtRgb)
	{
		RgbSpectrum rgb = UtilPbrt::RgbSpectrumFromPbrtVec3f(pbrtRgb);
		return make_shared<ConstantTexture<Spectrum>>(rgb);
	}

	static shared_ptr<Texture<Spectrum>> ParseTexture(const pbrt::Texture::SP & texture)
	{
		if (pbrt::ConstantTexture::SP pbrtConstantTexture = std::dynamic_pointer_cast<pbrt::ConstantTexture>(texture))
		{
			RgbSpectrum rgb = UtilPbrt::RgbSpectrumFromPbrtVec3f(pbrtConstantTexture->value);
			return make_shared<ConstantTexture<Spectrum>>(rgb);
		}
		else
		{
			std::cout << texture->toString() << std::endl;
		}

		return nullptr;
	}

	static shared_ptr<Bsdf> ParseBsdf(const pbrt::Material::SP & material)
	{
		if (pbrt::MatteMaterial::SP pbrtMatteMaterial = std::dynamic_pointer_cast<pbrt::MatteMaterial>(material))
		{
			shared_ptr<Texture<Spectrum>> texture;
			if (pbrtMatteMaterial->map_kd == nullptr)
			{
				texture = ParseTexture(pbrtMatteMaterial->kd);
			}
			else
			{
				texture = ParseTexture(pbrtMatteMaterial->map_kd);
			}
			shared_ptr<Bsdf> bsdf = make_shared<LambertDiffuseBsdf>(texture);
			return make_shared<TwoSidedBsdf>(bsdf);
		}
		else
		{
			shared_ptr<Bsdf> a = make_shared<LambertDiffuseBsdf>(make_shared<ConstantTexture<Spectrum>>(Spectrum(0.5)));
			return make_shared<TwoSidedBsdf>(a);
		}

		return nullptr;
	}

	static shared_ptr<Light> ParseLightSource(const pbrt::AreaLight::SP & pbrtLight)
	{
		if (pbrt::DiffuseAreaLightRGB::SP light = std::dynamic_pointer_cast<pbrt::DiffuseAreaLightRGB>(pbrtLight))
		{
			shared_ptr<const ConstantTexture<Spectrum>> emission = make_shared<const ConstantTexture<Spectrum>>(SpectrumConvertUtil::SpectrumFromRgb(light->L.x, light->L.y, light->L.z));
			return make_shared<DiffuseAreaLight>(emission);
		}
		else if (pbrt::DiffuseAreaLightBB::SP light = std::dynamic_pointer_cast<pbrt::DiffuseAreaLightBB>(pbrtLight))
		{
			shared_ptr<const ConstantTexture<Spectrum>> emission = make_shared<const ConstantTexture<Spectrum>>(SpectrumConvertUtil::SpectrumFromRgb(1.0, 1.0, 1.0));
			return make_shared<DiffuseAreaLight>(emission);
			//throw std::runtime_error("not support blackbody");
		}

		return nullptr;
	}

	struct ParsedGeometries
	{
		std::vector<shared_ptr<const Geometry>> mGeometries;
		std::vector<shared_ptr<const Geometry>> mLightGeometries;
		shared_ptr<EmbreeAccel> mAccel;
	};

	static ParsedGeometries ParseGeometries(const pbrt::Object::SP & object)
	{
		// create map for caching texture and material
		std::map<pbrt::MatteMaterial::SP, shared_ptr<Bsdf>> bsdfMap;
		std::map<pbrt::Texture::SP, shared_ptr<Texture<Spectrum>>> textureMap;

		// create result
		ParsedGeometries result;
		result.mAccel = make_shared<EmbreeAccel>();

		shared_ptr<EmbreeAccel> accel = result.mAccel;
		for (const pbrt::Shape::SP & shape : object->shapes)
		{
			std::vector<shared_ptr<const Geometry>> * geometries = (shape->areaLight == nullptr) ?
				&result.mGeometries : &result.mLightGeometries;

			shared_ptr<Geometry> geometry;

			// deal with geometry
			if (pbrt::TriangleMesh::SP mesh = std::dynamic_pointer_cast<pbrt::TriangleMesh>(shape))
			{
				// parse triangle mesh
				shared_ptr<TriangleMeshGeometry> triMesh = make_shared<TriangleMeshGeometry>();
				triMesh->mPositions.reserve(mesh->vertex.size());
				triMesh->mNormals.reserve(mesh->vertex.size());
				triMesh->mTextureCoords.reserve(mesh->vertex.size());
				for (int i = 0;i < mesh->vertex.size();i++)
				{
					const pbrt::vec3f & position = mesh->vertex[i];
					const pbrt::vec3f & normal = mesh->normal[i];
					triMesh->mPositions.emplace_back(position.x, position.y, position.z);
					if (i < mesh->normal.size()) triMesh->mNormals.emplace_back(Math::Normalize(Vec3(normal.x, normal.y, normal.z)));
					triMesh->mTextureCoords.emplace_back(0.0, 0.0);
				}
				for (int i = 0; i < mesh->index.size(); i++)
				{
					const pbrt::vec3i & face = mesh->index[i];
					triMesh->mTriangles.emplace_back(face.x, face.y, face.z);
				}

				triMesh->mBsdf = ParseBsdf(mesh->material);

				// add into geometry and accel
				geometry = triMesh;
				geometries->emplace_back(triMesh);
				accel->addGeometry(triMesh);
			}
			else
			{
				throw std::runtime_error("unknown shape type");
			}

			// deal with lightsource
			if (geometries == &result.mLightGeometries)
			{
				geometry->mAreaLight = ParseLightSource(shape->areaLight);
				geometry->precomputeForSampling();

				result.mLightGeometries.push_back(geometry);
			}

		}
		result.mAccel->commit();
		return result;
	}

	static shared_ptr<const Camera> ParseCamera(const pbrt::Camera::SP & camera, const Ivec2 & resolution)
	{
		return ThinlensCamera::FromPbrt(*camera, resolution);
	}

	static void Parse(const filesystem::path & pbrtFilepath)
	{
		pbrt::Scene::SP pbrtScene;
		if (pbrtFilepath.extension() == ".pbrt")
		{
			pbrtScene = pbrt::importPBRT(pbrtFilepath.string());
		}
		else if (pbrtFilepath.extension() == ".pbf")
		{
			pbrtScene = pbrt::Scene::loadFrom(pbrtFilepath.string());
		}
		else
		{
			throw std::runtime_error("unknown file extension");
		}
		
		// load camera
		shared_ptr<const Camera> camera = ParseCamera(pbrtScene->cameras[0], UtilPbrt::Ivec2FromPbrtVec2i(pbrtScene->film->resolution));

		// load geometries
		ParsedGeometries parsedGeometries = ParseGeometries(pbrtScene->world);

		// make scene from lightsource and meshes
		shared_ptr<const Scene> scene = make_shared<const Scene>(parsedGeometries.mAccel, parsedGeometries.mGeometries, parsedGeometries.mLightGeometries);

		shared_ptr<Sampler> sampler = make_shared<RandomSampler>(0);

		//shared_ptr<SingleCameraRenderer> renderer = make_shared<VisualizeGeometryNormal>(camera, scene, sampler, 10);
		//shared_ptr<SingleCameraRenderer> renderer = make_shared<PrimitivePathRenderer>(camera, scene, sampler, 2, 10000);
		shared_ptr<SingleCameraRenderer> renderer = make_shared<PathNeeRenderer>(camera, scene, sampler, 2, 10000);
		//shared_ptr<SingleCameraRenderer> renderer = make_shared<BidirPathRenderer>(camera, scene, sampler, 2, 10000);

		Viewer::Begin(camera->mFilm);
		renderer->render();
		Viewer::End(camera->mFilm);

		FimageIo::Save(camera->mFilm->getImage(), "test.pfm");
	}
};
