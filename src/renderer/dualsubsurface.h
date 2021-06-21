#pragma once

#include "common/wurst.h"

#include "common/ltmath.h"
#include "common/camera.h"
#include "common/renderer.h"
#include "common/scene.h"
#include "common/sampler.h"
#include "common/parallel.h"
#include "common/viewer.h"
#include "common/trimeshvertexiterator.h"
#include "sceneio/meshio/assimpio.h"

struct DualSubsurfaceRenderer : public SingleCameraRenderer
{
    static shared_ptr<DualSubsurfaceRenderer> FromSimple2Json(shared_ptr<const Camera> & camera,
                                                              shared_ptr<const Scene> & scene,
                                                              shared_ptr<Sampler> & sampler,
                                                              const json & json,
                                                              const filesystem::path & path)
    {
        return make_shared<DualSubsurfaceRenderer>(camera, scene, sampler);
    }

    DualSubsurfaceRenderer(shared_ptr<const Camera> & camera,
            shared_ptr<const Scene> & scene,
            shared_ptr<Sampler> & sampler):
            SingleCameraRenderer(camera, scene, sampler),
            mTriMeshVisitor(scene)
    {
    }

    void render() override
    {
        // illumination point sampling routine
        std::vector<Vec3> debugPositions;
        std::vector<SimpleVertex> points;
        std::vector<int> numPoints;

        // firstly, compute how many vertices we need
        mTriMeshVisitor.forEachMesh([&](int iMesh, shared_ptr<const TriangleMeshGeometry> & triMesh)
        {
            // precompute for sampling
            const_pointer_cast<TriangleMeshGeometry>(triMesh)->precomputeForSampling();
            // then for each mesh, samples some number of illumination points
            for (int iPoint = 0;iPoint < mNumPoints;iPoint++)
            {
                SimpleVertex sv = triMesh->sample(nullptr, mSampler->get2d());
                points.emplace_back(sv);
                debugPositions.push_back(sv.mPosition);
            }
        });

        // write the sampled point out just in case I made mistakes
        AssimpInOut::SaveObj("debug.obj", debugPositions);

        // illuminate the illumination point using monte carlo estimates
        for (int iPoint = 0; iPoint < points.size(); iPoint++)
        {
            // compute dipole approximation
            
        }
    }

    TriangleMeshVerticesVisitor mTriMeshVisitor;
    int mNumPoints = 1000;
};
