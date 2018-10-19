#define _CRT_SECURE_NO_DEPRECATE
#include "floatimage.h"

#include <malloc.h>
#include <cassert>
#include <memory>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

#include "common/floatimage/rgbe.h"

#include "ext/stb/stb_image.h"
#include "ext/stb/stb_image_write.h"

Fimage Fimage::LoadPFM(const std::string & filepath)
{
	std::ifstream is(filepath, std::ios::binary | std::ios::in);
    assert(is.is_open());
	size_t height;
	size_t width;
	std::string pfm;
	std::string minusOne;

	// read pfm header
	is >> pfm;
	is >> width >> height;
	is >> minusOne;
	assert(pfm == "PF");
	assert(width > 0);
	assert(height > 0);
	char c = is.get();
	if (c == '\r') c = is.get();

	// compute numpixels
	size_t numPixels = width * height;
	Fimage result(width, height);
	assert(false);


	// read data flip row
	for (size_t i = 0;i < height;i++)
	{
		std::vector<float> row(width * 3);
		is.read((char*)(&row[0]), sizeof(float) * width * 3);
	}

	is.close();

	return result;
}

void Fimage::SavePFM(const Fimage & fimage, const std::string & filepath)
{
	size_t width = fimage.mSize[0];
	size_t height = fimage.mSize[1];
	size_t numPixels = height * width;

	// write PFM header
	std::ofstream os(filepath, std::ios::binary);
	assert(os.is_open());
	os << "PF" << std::endl;
	os << width << " " << height << std::endl;
	os << "-1" << std::endl;

	// write data flip row
	for (size_t y = height;y >= 1;y--)
		for (size_t x = 0; x < width;x++)
			for (size_t i = 0; i < 3; i++)
			{
				const Vec3 v = fimage.colorAt(x, y - 1);
				float vf = static_cast<float>(v[i]);
				os.write((char*)(&vf), sizeof(float));
			}

	os.close();
}

Fimage Fimage::LoadHDR(const std::string & filepath)
{
	// open file
	FILE *fp = nullptr;
	fp = fopen(filepath.c_str(), "rb");
	assert(fp != nullptr);
	
	int width = 0;
	int height = 0;
	rgbe_header_info headerInfo;
	RGBE_ReadHeader(fp, &width, &height, &headerInfo);
	int numPixels = width * height;

	// read into hdrData
	Fimage result(width, height);
	std::vector<float> floats(width * height);
	RGBE_ReadPixels_RLE(fp, &(floats[0]), width, height);
	fclose(fp);

	result.setFloats(floats);

	return result;
}

Fimage Fimage::Load(const std::string & filepath)
{
	size_t i = filepath.find_last_of('.');
	assert(i > 0 && i < filepath.length() - 1);
	std::string extension = filepath.substr(i + 1);
	if (extension == "pfm")
		return LoadPFM(filepath);
	else if (extension == "hdr")
		return LoadHDR(filepath);
	else
		throw std::exception("unsupported file format");
}

void Fimage::SaveHDR(const Fimage & image, const std::string & filepath)
{
	FILE *fp = fopen(filepath.c_str(), "wb");
	assert(fp != nullptr);
    std::vector<float> floats = image.getFloats();
	RGBE_WriteHeader(fp, static_cast<int>(image.mSize[0]), static_cast<int>(image.mSize[1]), NULL);
	RGBE_WritePixels_RLE(fp, &floats[0], static_cast<int>(image.mSize[0]), static_cast<int>(image.mSize[1]));
	fclose(fp);
}

void Fimage::SavePNG(const Fimage & fimage, const std::string & filepath)
{
	char * data = new char[fimage.mSize[0] * fimage.mSize[1] * 3];
	for (size_t row = 0;row < fimage.mSize[1];row++)
		for (size_t col = 0;col < fimage.mSize[0];col++)
			for (size_t i = 0;i < 3;i++)
			{
				double p = std::pow(fimage.colorAt(col, row)[i], double(1 / 2.2));
				p = std::min(p * 255.9999, 255.0);
				data[(col + row * fimage.mSize[0]) * 3 + i] = static_cast<char>(p);
			}
	stbi_write_png(filepath.c_str(),
				   static_cast<int>(fimage.mSize[0]),
				   static_cast<int>(fimage.mSize[1]),
				   3,
				   data,
				   static_cast<int>(fimage.mSize[0] * 3));
    delete data;
}

void Fimage::Save(const Fimage & fimage, const std::string & filepath)
{
	size_t i = filepath.find_last_of('.');
	assert(i > 0 && i < filepath.length() - 1);
	std::string extension = filepath.substr(i + 1);
	if (extension == "pfm")
		SavePFM(fimage, filepath);
	else if (extension == "hdr")
		SaveHDR(fimage, filepath);
	else if (extension == "png")
		SavePNG(fimage, filepath);
	else
		throw std::exception("unsupported file format");
}

Fimage Fimage::Pow(const Fimage & image, const float exponent)
{
	size_t numRows = image.mSize[1];
	size_t numCols = image.mSize[0];
	Fimage result(image.mSize);

	for (size_t row = 0;row < numRows;row++)
		for (size_t col = 0;col < numCols;col++)
		{
			result.mData[col + row * numCols] = Vec3::Pow(image.mData[col + row * numCols], exponent);
		}
	return result;
}

Fimage Fimage::ResizeHalf(const Fimage & image)
{
	assert(image.mSize[0] % 2 == 0 && image.mSize[1] % 2 == 0);
	Fimage result(image.mSize[0] / 2, image.mSize[1] / 2);
	result.forEachPixel([&] (const size_t x, const size_t y, Vec3 * color) {
		Vec3 tl = image.colorAt(x * 2, y * 2);
		Vec3 tr = image.colorAt(x * 2 + 1, y * 2);
		Vec3 bl = image.colorAt(x * 2, y * 2 + 1);
		Vec3 br = image.colorAt(x * 2 + 1, y * 2 + 1);
		*color = (tl + tr + bl + br) / 4.0f;
	});
	return result;
}

Fimage::Fimage(const size_t sizeX, const size_t sizeY, const std::vector<Vec3>& data):
	mSize(Uvec2(sizeX, sizeY)),
	mData(data),
	mWrapS(Fimage::WrapMode::Clamp),
	mWrapT(Fimage::WrapMode::Clamp)
{
	assert(data.size() == sizeX * sizeY);
}

Fimage::Fimage(const size_t sizeX, const size_t sizeY) :
	mSize(Uvec2(sizeX, sizeY)),
	mData(sizeX * sizeY),
	mWrapS(Fimage::WrapMode::Clamp),
	mWrapT(Fimage::WrapMode::Clamp)
{
}

Fimage::Fimage(const Uvec2 & size) :
	Fimage(size[0], size[1])
{
}

void Fimage::setZero()
{
	std::memset(&mData[0][0], 0, sizeof(Vec3) * mSize[0] * mSize[1]);
}

Vec3 Fimage::evalTexel(int32_t x, int32_t y) const
{
	const Uvec2 & size = mSize;
	switch (mWrapS)
	{
		case Fimage::Repeat:
			x = Math::PositiveMod(x, static_cast<int32_t>(size[0]));
			break;
		case Fimage::Clamp:
			x = Math::Clamp(x, 0, static_cast<int32_t>(size[0] - 1));
			break;
		case Fimage::Mirror:
			assert(false); /// TODO:: Implement
			break;
		default:
			assert(false);
			break;
	}

	switch (mWrapT)
	{
		case Fimage::Repeat:
			y = Math::PositiveMod(y, static_cast<int32_t>(size[1]));
			break;
		case Fimage::Clamp:
			y = Math::Clamp(y, 0, static_cast<int32_t>(size[1] - 1));
			break;
		case Fimage::Mirror:
			assert(false); /// TODO:: Implement
			break;
		default:
			assert(false);
			break;
	}

	return this->colorAt(x, y);
}

Vec3 Fimage::evalNearest(const Vec2 & uv) const
{
	// taken from mitsuba mipmap.h
	const Vec2 &size = mSize;
	return evalTexel(Math::FloorToInt(uv[0] * size[0]), Math::FloorToInt(uv[1] * size[1]));
}

Vec3 Fimage::evalBilinear(const Vec2 & uv) const
{
	// taken from mitsuba mipmap.h
	const Vec2 &size = mSize;
	Vec2 scaledUv = uv * size - Vec2(0.5);

	int32_t xPos = Math::FloorToInt(scaledUv[0]), yPos = Math::FloorToInt(scaledUv[1]);
	double dx1 = scaledUv[0] - xPos, dy1 = scaledUv[1] - yPos;
	double dx2 = 1.0 - dx1, dy2 = 1.0 - dy1;

	return evalTexel(xPos, yPos) * dx2 * dy2
		+ evalTexel(xPos, yPos + 1) * dx2 * dy1
		+ evalTexel(xPos + 1, yPos) * dx1 * dy2
		+ evalTexel(xPos + 1, yPos + 1) * dx1 * dy1;
}

void Fimage::forEachPixel(const std::function<void(const size_t x, const size_t y, Vec3*color)>& func)
{
	for (size_t y = 0;y < mSize[1];y++)
	{
		const size_t xOffset = y * mSize[0];
		for (size_t x = 0;x < mSize[0];x++)
		{
			func(x, y, &mData[xOffset + x]);
		}
	}
}

void Fimage::forEachPixel(const std::function<void(const size_t x, const size_t y, const Vec3&color)>& func) const
{
	for (size_t y = 0; y < mSize[1]; y++)
	{
		const size_t xOffset = y * mSize[0];
		for (size_t x = 0; x < mSize[0]; x++)
		{
			func(x, y, mData[xOffset + x]);
		}
	}
}

void Fimage::forEachPixel(const std::function<void(const Vec2 & uv, Vec3 * color)> & func)
{
	Vec2 invResolution = 1.0f / Vec2(mSize);
	std::cout << invResolution << std::endl;
	for (size_t y = 0;y < mSize[1];y++)
	{
		const size_t xOffset = y * mSize[0];
		for (size_t x = 0;x < mSize[0];x++)
		{
			Vec2 uv = Vec2(x + 0.5f, y + 0.5f) * invResolution;
			func(uv, &mData[xOffset + x]);
		}
	}
}

void Fimage::forEachPixel(const std::function<void(const Vec2 & uv, const Vec3 & color)> & func) const
{
	Vec2 invResolution = 1.0f / Vec2(mSize);
	for (size_t y = 0;y < mSize[1];y++)
	{
		const size_t xOffset = y * mSize[0];
		for (size_t x = 0;x < mSize[0];x++)
		{
			Vec2 uv = Vec2(x + 0.5f, y + 0.5f) * invResolution;
			func(uv, mData[xOffset + x]);
		}
	}
}

void Fimage::setFloats(std::vector<float> floats)
{
	forEachPixel([&](const size_t x, const size_t y, Vec3 * color) {
		for (size_t i = 0; i < 3; i++)
			(*color)[i] = floats[(y * mSize[0] + x) * 3 + i];
	});
}

const std::vector<float> Fimage::getFloats() const
{
	std::vector<float> result;
	forEachPixel([&](const size_t x, const size_t y, const Vec3 & color) {
		for (size_t i = 0; i < 3; i++)
			result.emplace_back(static_cast<float>(color[i]));
	});
    return result;
}

Fimage Fimage::operate(const Fimage & image, const std::function<Vec3(const Vec3&in, const Vec3&out)>& func) const
{
	assert(image.mSize == this->mSize);
	Fimage result(this->mSize);
	for (size_t y = 0;y < this->mSize[1];y++)
		for (size_t x = 0;x < this->mSize[0];x++)
			result.colorAt(x, y) = func(this->colorAt(x, y), image.colorAt(x, y));
	return result;
}

Fimage Fimage::operator+(const Fimage & image) const
{
	return operate(image, [](const Vec3 & in, const Vec3 & out) { return in + out; });
}

Fimage Fimage::operator*(const Fimage & image) const
{
	return operate(image, [](const Vec3 & in, const Vec3 & out) { return in * out; });
}
