#pragma once

#include "common/math/math.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <functional>

class Fimage
{
public:
	enum WrapMode
	{
		Repeat,
		Clamp,
		Mirror,
	};

	static Fimage LoadPFM(const std::string &filepath);
	static Fimage LoadHDR(const std::string &filepath);
	static Fimage Load(const std::string &filepath);
	static void SavePFM(const Fimage &fimage, const std::string &filepath);
	static void SaveHDR(const Fimage &image, const std::string &filepath);
	static void SavePNG(const Fimage &fimage, const std::string &filepath);
	static void Save(const Fimage &fimage, const std::string &filepath);
	static Fimage Pow(const Fimage &image, const float value);
	static Fimage ResizeHalf(const Fimage & image);

	Fimage() : mSize(Uvec2(0, 0)) {}
	Fimage(const size_t sizeX, const size_t sizeY, const std::vector<Vec3> &data);
	Fimage(const size_t sizeX, const size_t sizeY);
	Fimage(const Uvec2 &size);

	void setZero();

	inline Vec3& colorAt(const size_t x, const size_t y) { assert(x < mSize[0]); assert(y < mSize[1]); return mData[y * mSize[0] + x]; }
	inline const Vec3& colorAt(const size_t x, const size_t y) const { assert(x < mSize[0]); assert(y < mSize[1]); return mData[y * mSize[0] + x]; }

	Vec3 evalTexel(int32_t x, int32_t y) const;
	Vec3 evalNearest(const Vec2 & uv) const;
	Vec3 evalBilinear(const Vec2 & uv) const;

	void forEachPixel(const std::function<void(const size_t x, const size_t y, Vec3 * color)> & func);
	void forEachPixel(const std::function<void(const size_t x, const size_t y, const Vec3 & color)> & func) const;
	void forEachPixel(const std::function<void(const Vec2 & uv, Vec3 * color)> & func);
	void forEachPixel(const std::function<void(const Vec2 & uv, const Vec3 & color)> & func) const;

	inline const size_t getNumPixels() const { return mSize[0] * mSize[1]; }

	void setFloats(std::vector<float> floats);
	const std::vector<float> getFloats() const;

	Fimage operate(const Fimage & image, const std::function<Vec3(const Vec3 & in, const Vec3 & out)> & func) const;
	Fimage operator+(const Fimage & image) const;
	Fimage operator*(const Fimage & image) const;

	WrapMode mWrapS = WrapMode::Repeat;
	WrapMode mWrapT = WrapMode::Repeat;
	std::vector<Vec3> mData;
	Uvec2 mSize;
};
