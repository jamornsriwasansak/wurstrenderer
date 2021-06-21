#pragma once

#include "common/wurst.h"

#include "common/scene.h"
#include "common/camera.h"
#include "common/viewer.h"

#include "accel/embreeaccel.h"
#include "bsdf/lambertdiffuse.h"
#include "bsdf/microfacet.h"
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

struct TungstenParser
{
	using BsdfDict = std::map<const std::string, shared_ptr<const Bsdf>>;

	static shared_ptr<const Bsdf> ParseBsdf(const json & bsdfJson, const filesystem::path & abspath)
	{
		const std::string & type = bsdfJson["type"];

		if (type == "lambert")
		{ 
			return make_shared<const LambertDiffuseBsdf>(ParseTexture(bsdfJson["albedo"], abspath));
		}

		if (type == "oren_nayar")
		{
			// TODO:: multiplied roughness by 1.0 / std::sqrt(2)
			return make_shared<const OrenNayarBsdf>(ParseTexture(bsdfJson["albedo"], abspath), /*1.0 / std::sqrt(2) */ ParseTexture(bsdfJson["roughness"], abspath));
		}

		if (type == "phong")
		{
			// TODO::
			//return make_shared<const PhongBsdf>(ParseTexture(bsdfJson["albedo"], abspath), ParseTexture(bsdfJson["exponent"], abspath));
		}

		if (type == "plastic")
		{
			// TODO::
		}

		if (type == "rough_plastic")
		{
			// TODO::
		}

		if (type == "rough_conductor")
		{
			// TODO::
			//return make_shared<const MicrofacetReflection>(ParseTexture(bsdfJson["roughness"], abspath), ParseTexture(bsdfJson["roughness"], abspath));
			//return make_shared<const MicrofacetReflectionBsdf>(make_shared<ConstantTexture>(Spectrum(1.0)), make_shared<ConstantTexture>(Spectrum(1.0)));
		}

		if (type == "smooth_coat")
		{
			//return make_shared<const MicrofacetReflectionBsdf>(make_shared<ConstantTexture>(Spectrum(0.1)), make_shared<ConstantTexture>(Spectrum(0.1)));
		}

		if (type == "null")
		{
			return make_shared<const LambertDiffuseBsdf>(ParseTexture(bsdfJson["albedo"], abspath));
		}

		std::cout << "encounter unknown bsdf : " << type << std::endl;

		return make_shared<const LambertDiffuseBsdf>(make_shared<ConstantTexture<Spectrum>>(Spectrum(0.5)));
		//return make_shared<const OrenNayarBsdf>(make_shared<ConstantTexture>(Spectrum(0.5)), make_shared<ConstantTexture>(Spectrum(1.0)));
	}

	static BsdfDict ParseBsdfs(const json & bsdfJsons, bool isTwoSided, const filesystem::path & abspath)
	{
		BsdfDict result;
		for (const auto & bsdfJson : bsdfJsons)
		{
			const std::string name = bsdfJson["name"];
			shared_ptr<const Bsdf> bsdf = ParseBsdf(bsdfJson, abspath);
			if (bsdf != nullptr)
			{
				bsdf = isTwoSided ? make_shared<const TwoSidedBsdf>(bsdf) : bsdf;
			}
			result[name] = bsdf;
		} 
		return result;
	}

	static shared_ptr<const Camera> ParseCamera(const json & cameraJson)
	{
		Ivec2 resolution = Ivec2FromJson(cameraJson["resolution"]);
		
		// assume type = thinlens
		Vec3 position = Vec3FromJson(cameraJson["transform"]["position"]);
		Vec3 lookAt = Vec3FromJson(cameraJson["transform"]["look_at"]);
		Vec3 up = Vec3FromJson(cameraJson["transform"]["up"]);

		double fovy = double(cameraJson["fov"]) * Math::Pi / 180.0 / double(resolution[0]) * double(resolution[1]);

		// TODO:: fix this random parameter
		return make_shared<const ThinlensCamera>(position, lookAt, up, fovy, resolution, 0.01, 6.5);
	}

	static shared_ptr<DiffuseAreaLight> ParseArealight(const json & emissionJson)
	{
		return make_shared<DiffuseAreaLight>(SpectrumConvertUtil::SpectrumsFromRgbs(Vec3FromJson(emissionJson)));
	}

	static shared_ptr<EnvmapLight> ParseEnvmapLight(const json & emissionJson, const filesystem::path & abspath)
	{
		if (emissionJson.type() == json::value_t::string)
		{
			const std::string envmapFilepath = emissionJson;
			shared_ptr<SphRep<Spectrum>> sphrep = make_shared<EquiRectRep<Spectrum>>(FimageIo::Load(nullptr, abspath.parent_path() / envmapFilepath), true);
			return make_shared<EnvmapLight>(sphrep);
		}
		else
		{
			shared_ptr<SphRep<Spectrum>> sphrep = make_shared<EquiRectRep<Spectrum>>(SpectrumConvertUtil::SpectrumsFromRgbs(Vec3FromJson(emissionJson)));
			return make_shared<EnvmapLight>(sphrep);
		}
	}

	static shared_ptr<Scene> ParseGeometry(const json & primitivesJson, const BsdfDict & bsdfDict, const bool isTwoSided, const filesystem::path & abspath)
	{
		shared_ptr<EmbreeAccel> accel = make_shared<EmbreeAccel>();
		shared_ptr<EnvmapLight> envmapLight = nullptr;
		std::vector<shared_ptr<const Geometry>> geometries;
		std::vector<shared_ptr<const Geometry>> lightGeometries;

		for (const auto & primitiveJson : primitivesJson)
		{
			// decide whether this geometry is surface or arealight
			const std::string & typeStr = primitiveJson["type"];
			shared_ptr<const Bsdf> bsdf = nullptr;
			shared_ptr<Light> arealight = nullptr;

			if (primitiveJson.has("emission"))
			{
				if (typeStr == "infinite_sphere")
				{
					if (envmapLight != nullptr) throw std::runtime_error("should be only 1 envmap");
					envmapLight = ParseEnvmapLight(primitiveJson["emission"], abspath);
				}
				else
				{
					arealight = ParseArealight(primitiveJson["emission"]);
				}
			}
			else
			{
				if (primitiveJson["bsdf"].type() == json::value_t::object)
				{
					bsdf = ParseBsdf(primitiveJson["bsdf"], abspath);
				}
				else
				{
					const std::string & bsdfStr = primitiveJson["bsdf"];
					if (bsdfDict.find(bsdfStr) != bsdfDict.end())
					{
						bsdf = bsdfDict.at(bsdfStr);
					}
					else
					{
						bsdf = make_shared<const LambertDiffuseBsdf>(make_shared<ConstantTexture<Spectrum>>(Spectrum(0.5)));
						bsdf = isTwoSided ? make_shared<const TwoSidedBsdf>(bsdf) : bsdf;
					}
				}
			}

			Matrix4 transform;
			if (primitiveJson.has("transform"))
			{
				transform = ParseTransform(primitiveJson["transform"]);
			}

			shared_ptr<Geometry> geometry;
			if (typeStr == "mesh")
			{
				const std::string & file = primitiveJson["file"];
				shared_ptr<TriangleMeshGeometry> tg = Wo3InOut::LoadFromWo3(abspath.parent_path() / file);
				geometry = tg;
				accel->addGeometry(tg);
			}
			else if (typeStr == "quad")
			{
				shared_ptr<PlaneGeometry> pg = make_shared<PlaneGeometry>(Vec3(0), Vec3(0, 0, 1), Vec3(1, 0, 0));
				geometry = pg;
				accel->addGeometry(pg);
			}
			else if (typeStr == "cube")
			{
				shared_ptr<CubeGeometry> cg = make_shared<CubeGeometry>();
				geometry = cg;
				accel->addGeometry(cg);
			}

			geometry->mBsdf = bsdf;
			if (arealight)
			{
				geometry->mAreaLight = arealight;
				geometry->transform(transform);
				geometry->precomputeForSampling();
				lightGeometries.emplace_back(geometry);
			}
			geometries.push_back(geometry);
		}
		accel->commit();
		return make_shared<Scene>(accel, geometries, lightGeometries, envmapLight);
	}

	static shared_ptr<Renderer> ParseRenderer(const json & integratorJson,
													const json & rendererJson,
													shared_ptr<const Camera> & camera,
													shared_ptr<const Scene> & scene,
													shared_ptr<Sampler> & sampler)
	{
		// read num bounces
		int maxBounces = integratorJson.getValue("max_bounces", 32);
		assert(maxBounces >= 0);
		int numVertices = maxBounces + 2;
		
		// read num spp
		int spp = rendererJson.getValue("spp", 64);
		assert(spp >= 0);
		int numSamples = spp;

		std::string integratorName = integratorJson.getValue("type", "primitive_path_tracer");
		if (integratorName == "primitive_path_tracer")
		{
			return make_shared<PrimitivePathRenderer>(camera, scene, sampler, numVertices, numSamples);
		}
		else if (integratorName == "path_tracer")
		{
			return make_shared<PathNeeRenderer>(camera, scene, sampler, numVertices, numSamples);
		}
		else if (integratorName == "light_tracer")
		{
			return make_shared<LightSplattingRenderer>(camera, scene, sampler, numVertices, numSamples);
		}
		else if (integratorName == "kelemen_mlt")
		{
			return make_shared<PssmltPrimitivePathRenderer>(camera, scene, sampler, numVertices, numSamples);
		}
		else if (integratorName == "photon_map")
		{
			return make_shared<PhotonMapRenderer>(camera, scene, sampler, numVertices, numSamples);
		}

        return make_shared<VisualizeGeometryNormal>(camera, scene, sampler, numSamples);
	}

	static shared_ptr<const Texture<Spectrum>> ParseTexture(const json & albedoJson, const filesystem::path & abspath)
	{
		if (albedoJson.type() == json::value_t::string)
		{
			std::string pathToImage = albedoJson;
			//return make_shared<const ImageTexture<Spectrum>>(abspath.parent_path() / pathToImage);
			// TODO::
			return nullptr;
		}

		if (!albedoJson.has("type"))
		{
			Spectrum albedo = SpectrumConvertUtil::SpectrumsFromRgbs(Vec3FromJson(albedoJson));
			return make_shared<const ConstantTexture<Spectrum>>(albedo);
		}

		const std::string & type = albedoJson["type"];
		if (type == "checker")
		{
			Spectrum color1 = SpectrumConvertUtil::SpectrumsFromRgbs(Vec3FromJson(albedoJson["on_color"]));
			Spectrum color2 = SpectrumConvertUtil::SpectrumsFromRgbs(Vec3FromJson(albedoJson["off_color"]));
			Uint resU = albedoJson["res_u"];
			Uint resV = albedoJson["res_v"];
			return make_shared<CheckerboardTexture<Spectrum>>(color1, color2, Uvec2(resU, resV));
		}

		return nullptr;
	}

	static Matrix4 ParseTransform(const json & transformJson)
	{
		Matrix4 transform;
		if (transformJson.has("scale")) { transform = Matrix4::Scale(Vec3FromJson(transformJson["scale"])) * transform; }
		if (transformJson.has("rotation")) { transform = Matrix4::RotateYXZ(Math::RadiansFromDegrees(Vec3FromJson(transformJson["rotation"]))) * transform; }
		if (transformJson.has("position")) { transform = Matrix4::Translate(Vec3FromJson(transformJson["position"])) * transform; }
		return transform;
	}

	static void Parse(const filesystem::path & filepath)
	{
		// read a file
		std::ifstream ifs(filepath);
		if (!ifs.is_open()) throw std::runtime_error((std::string("file not found") + filepath.string()).c_str());
		json sceneJson;
		ifs >> sceneJson;

		filesystem::path abspath = Util::File::MakeAbsolutePath(filepath);

		// read all bsdfs
		const BsdfDict bsdfDict = (sceneJson.has("bsdfs")) ? ParseBsdfs(sceneJson["bsdfs"], true, abspath) : BsdfDict();

		// get camera
		shared_ptr<const Camera> camera = ParseCamera(sceneJson["camera"]);

		// get scene file
		shared_ptr<const Scene> scene = ParseGeometry(sceneJson["primitives"], bsdfDict, true, abspath);

		// TODO:: get sampler
		shared_ptr<Sampler> sampler = make_shared<RandomSampler>(0);

		// get renderer
		shared_ptr<Renderer> renderer = ParseRenderer(sceneJson["integrator"], sceneJson["renderer"], camera, scene, sampler);

		// create viewer and run
		Viewer::Begin(camera->mFilm);

		// render!
		renderer->render();

		FimageIo::Save(camera->mFilm->getImage(), abspath.parent_path() / "wip.pfm");

		// close viewer
		Viewer::End(camera->mFilm);
	}
};