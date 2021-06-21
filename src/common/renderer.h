#pragma once

#include "common/wurst.h"

struct Renderer
{
	Renderer(shared_ptr<const Scene> & scene, shared_ptr<Sampler> & sampler):
		mScene(scene),
		mSampler(sampler)
	{
	}

	virtual void render() = 0;

	shared_ptr<Sampler> mSampler = nullptr;
	shared_ptr<const Scene> mScene = nullptr;
};

struct Splat
{
	Splat() {}

	Splat(const Vec2 & raster, const Spectrum & contrib) : mRaster(raster), mContrib(contrib) {}

	Vec2 mRaster;
	Spectrum mContrib;
};

struct SingleCameraRenderer : public Renderer
{
	SingleCameraRenderer(shared_ptr<const Camera> camera, shared_ptr<const Scene> scene, shared_ptr<Sampler> sampler):
		mCamera(camera),
		Renderer(scene, sampler)
	{
	}

	virtual void render() override = 0;

	shared_ptr<const Camera> mCamera = nullptr;
};

struct SingleCameraRenderer2 : public Renderer
{
	SingleCameraRenderer2(shared_ptr<const Camera2> camera, shared_ptr<const Scene> scene, shared_ptr<Sampler> sampler):
		mCamera2(camera),
		Renderer(scene, sampler)
	{
	}

	virtual void render() override = 0;

	shared_ptr<const Camera2> mCamera2 = nullptr;
};
