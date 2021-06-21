#pragma once

#include "common/wurst.h"
#include "sceneio/meshio/assimpio.h"

#include "common/scene.h"
#include "common/viewer.h"
#include "camera2/circlegrid.h"
#include "camera/thinlens.h"
#include "light/envmap.h"
#include "medium/homogeneous.h"
#include "renderer/pathnextevent.h"
#include "renderer/lightsplatting.h"
#include "renderer/photonmap.h"
#include "renderer/primitivepath.h"
#include "renderer/pathnextevent.h"
#include "renderer/pssmltprimitivepath.h"
#include "renderer/visualizegeomnormal.h"
#include "renderer/diffuseprt.h"
#include "renderer/glossyprt2002.h"
#include "renderer/dualsubsurface.h"
#include "renderer2/primitivepath2.h"
#include "sphericalrepresentation/equirectangular.h"

struct Simple2Parser
{
    static shared_ptr<const Scene> ParseScene(const json & meshesJson, const shared_ptr<EnvmapLight> & envmap, const filesystem::path & rootJsonFullPath)
    {
        // load meshes and feed into accel
        shared_ptr<EmbreeAccel> accel = make_shared<EmbreeAccel>();
        std::vector<shared_ptr<const Geometry>> geometries;
        std::vector<shared_ptr<const Geometry>> lightGeometries;
        for (const json & json : meshesJson)
        {
            const filesystem::path objPath = rootJsonFullPath.parent_path() / Util::Json::GetValue<std::string>(json, "path", "");
            const bool doSrgbConversionTheMitsubaWay = Util::Json::GetValue<bool>(json, "mitsuba_srgb", false);
            const bool doUsePhongTessellation = Util::Json::GetValue<bool>(json, "phong_tessellation", false);
            AssimpInOut::Result obj = AssimpInOut::LoadFromObj(objPath, doSrgbConversionTheMitsubaWay);
            for (const shared_ptr<TriangleMeshGeometry> & p : obj.mGeometries)
            {
                if (doUsePhongTessellation)
                {
                    shared_ptr<PhongMeshGeometry> q = make_shared<PhongMeshGeometry>();
                    std::swap(q->mPositions, p->mPositions);
                    std::swap(q->mNormals, p->mNormals);
                    std::swap(q->mTextureCoords, p->mTextureCoords);
                    std::swap(q->mTriangles, p->mTriangles);
                    q->mAlpha = p->mAlpha;
                    q->mBsdf = p->mBsdf;
                    q->mMediumBackface = p->mMediumBackface;
                    q->mMediumFrontface = p->mMediumFrontface;
                    accel->addUserGeometry(q);
                }
                else
                {
                    accel->addGeometry(p);
                }
            }
            // note: we don't have sampling routine in phong mesh tessellation. we thus use triangle geometry as a lightsource for now.
            for (const shared_ptr<TriangleMeshGeometry> & p : obj.mLightGeometries) accel->addGeometry(p);
            geometries.insert(geometries.end(), obj.mGeometries.begin(), obj.mGeometries.end());
            lightGeometries.insert(lightGeometries.end(), obj.mLightGeometries.begin(), obj.mLightGeometries.end());
        }
        accel->commit();
        // make scene from lightsource and meshes
        shared_ptr<const Scene> scene = make_shared<const Scene>(accel, geometries, lightGeometries, envmap);
        return scene;
    }

    static shared_ptr<const Medium> ParseMedium(const json & mediumJson)
    {
        if (mediumJson["type"] == "homogeneous")
        {
            return HomogeniusMedium::FromSimple2Json(mediumJson);
        }
        else if (mediumJson["type"] == "none")
        {
            return nullptr;
        }

        throw std::runtime_error("unknown medium");
    }

    static void Parse(const filesystem::path & rootJsonPath)
    {
        filesystem::path rootJsonFullPath = Util::File::MakeAbsolutePath(rootJsonPath);

        // read a file
        std::ifstream ifs(rootJsonFullPath, std::ifstream::in);
        if (!ifs.is_open()) throw std::runtime_error((std::string("file not found") + rootJsonFullPath.string()).c_str());
        json sceneJson;
        ifs >> sceneJson;
        Parse(sceneJson, rootJsonFullPath);
    }

    static void Parse(const json & sceneJson)
    {
        filesystem::path path = Util::File::MakeAbsolutePath(filesystem::current_path()) / "<none>.json";
        Parse(sceneJson, path);
    }

    static void Parse(const json & sceneJson, const filesystem::path & rootJsonPath)
    {
        filesystem::path rootJsonFullPath = Util::File::MakeAbsolutePath(rootJsonPath);

        // check comment and if use viewer
        const std::string comment = Util::Json::GetValue<std::string>(sceneJson, "comment", "");
        const bool useViewer = Util::Json::GetValue<bool>(sceneJson, "viewer", false);

        // load envmap
        shared_ptr<EnvmapLight> envmapLight = nullptr;
        {
            std::string envmapPath = Util::Json::GetValue<std::string>(sceneJson, "envmap", "");
            if (envmapPath != "")
            {
                Fimage<Spectrum> envmapImage = FimageIo::Load(rootJsonFullPath.parent_path() / envmapPath);
                shared_ptr<const SphRep<Spectrum>> envmapRep = make_shared<const EquiRectRep<Spectrum>>(envmapImage, true);
                envmapLight = make_shared<EnvmapLight>(envmapRep);
            }
        }

        // create sampler
        shared_ptr<Sampler> sampler = make_shared<RandomSampler>(0);
        shared_ptr<const Scene> scene = (sceneJson.find("meshes") == sceneJson.end()) ? nullptr : ParseScene(sceneJson["meshes"], envmapLight, rootJsonFullPath);

        // render 2d
        if (sceneJson.count("render2d") > 0)
        {
            shared_ptr<const Camera2> camera2 = CircleGridCamera2::FromSimple2ParseJson(sceneJson["camera2d"]);
            json rendersJson2d = sceneJson["render2d"];
            for (int iRender = 0; iRender < rendersJson2d.size(); iRender++)
            {
                json rendererJson = rendersJson2d[iRender];

                // create a renderer
                shared_ptr<Renderer> renderer = nullptr;
                if (rendererJson["integrator"] == "primitive_path_tracer2d")
                {
                    renderer = PrimitivePathRenderer2::FromSimple2Json(camera2, scene, sampler, rendererJson, rootJsonFullPath);
                }
                assert(renderer != nullptr);

                // read output filename
                const std::string outputFile = Util::Json::GetValue<std::string>(rendererJson, "output", "");

                // render and view if needed
                camera2->mFilm->mName = comment;
                if (useViewer) Viewer::Begin(camera2->mFilm);
                renderer->render();
                if (useViewer) Viewer::End(camera2->mFilm);

                // save
                if (outputFile != "") FimageIo::Save(camera2->mFilm->getImage(), rootJsonFullPath.parent_path() / outputFile);
            }
        }

        // render 3d
        if (sceneJson.count("render") > 0)
        {
            // camera
            shared_ptr<const Camera> camera = ThinlensCamera::FromSimple2ParseJson(sceneJson["camera"]);

            if (sceneJson.count("medium"))
            {
                shared_ptr<const Medium> main_medium = ParseMedium(sceneJson["medium"]);
                camera->mMedium = main_medium;
            }

            json rendersJson = sceneJson["render"];
            for (int iRender = 0; iRender < rendersJson.size(); iRender++)
            {
                json rendererJson = rendersJson[iRender];

                // create a renderer
                shared_ptr<Renderer> renderer = nullptr;
                if (rendererJson["integrator"] == "bidirectional_path_tracer")
                {
                    renderer = BidirPathRenderer::FromSimple2Json(camera, scene, sampler, rendererJson, rootJsonFullPath);
                }
                else if (rendererJson["integrator"] == "photon_map")
                {
                    renderer = PhotonMapRenderer::FromSimple2Json(camera, scene, sampler, rendererJson, rootJsonFullPath);
                }
                else if (rendererJson["integrator"] == "primitive_path_tracer")
                {
                    renderer = PrimitivePathRenderer::FromSimple2Json(camera, scene, sampler, rendererJson, rootJsonFullPath);
                }
                else if (rendererJson["integrator"] == "light_splatting")
                {
                    renderer = LightSplattingRenderer::FromSimple2Json(camera, scene, sampler, rendererJson, rootJsonFullPath);
                }
                else if (rendererJson["integrator"] == "nextevent_path_tracer")
                {
                    renderer = PathNeeRenderer::FromSimple2Json(camera, scene, sampler, rendererJson, rootJsonFullPath);
                }
                else if (rendererJson["integrator"] == "visualize")
                {
                    renderer = VisualizeGeometryNormal::FromSimple2Json(camera, scene, sampler, rendererJson, rootJsonFullPath);
                }
                else if (rendererJson["integrator"] == "diffuse_prt")
                {
                    renderer = DiffusePrt::FromSimple2Json(camera, scene, sampler, rendererJson, rootJsonFullPath);
                }
                else if (rendererJson["integrator"] == "glossy_prt2002")
                {
                    renderer = GlossyPrt2002::FromSimple2Json(camera, scene, sampler, rendererJson, rootJsonFullPath);
                }
                assert(renderer != nullptr);

                // read output filename
                const std::string outputFile = Util::Json::GetValue<std::string>(rendererJson, "output", "");

                // render and view if needed
                camera->mFilm->mName = comment;
                if (useViewer) Viewer::Begin(camera->mFilm);
                renderer->render();
                if (useViewer) Viewer::End(camera->mFilm);

                // save
                if (outputFile != "") FimageIo::Save(camera->mFilm->getImage(), rootJsonFullPath.parent_path() / outputFile);
            }
        }
    }
};
