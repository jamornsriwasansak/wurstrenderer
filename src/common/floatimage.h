#pragma once

#include "common/wurst.h"

#include "ext/stb/stb_image.h"
#include "ext/stb/stb_image_write.h"
#include "ext/tinyexr/tinyexr.h"

template <typename Type>
struct Fimage
{
	enum WrapMode
	{
		Repeat,
		Clamp,
		Mirror,
	};

	static Fimage FlipY(const Fimage & v)
	{
		Fimage result(v.mSize);
		for (int y = 0; y < v.mSize[1]; y++)
			for (int x = 0; x < v.mSize[0]; x++)
				result.at(x, y) = v.at(x, v.mSize[1] - y - 1);
		return result;
	}

	template <typename Scalar>
	Fimage(const std::vector<Scalar> & image, const int sizeX, const int sizeY): mData(sizeX * sizeY), mSize(sizeX, sizeY)
	{
		assert(image.size() == sizeX * sizeY);
		for (int y = 0; y < sizeY; y++)
			for (int x = 0; x < sizeX; x++)
				at(x, y) = Type(image[y * sizeX + x]);
	}

	Fimage(): mSize(Ivec2(0, 0))
	{
	}

	Fimage(const int sizeX, const int sizeY, const std::vector<Type> & data):
		mSize(Ivec2(sizeX, sizeY)),
		mData(data),
		mWrapS(Fimage::WrapMode::Clamp),
		mWrapT(Fimage::WrapMode::Clamp)
	{
		assert(data.size() == sizeX * sizeY);
	}

	Fimage(const int sizeX, const int sizeY):
		mSize(Ivec2(sizeX, sizeY)),
		mData(sizeX * sizeY),
		mWrapS(Fimage::WrapMode::Clamp),
		mWrapT(Fimage::WrapMode::Clamp)
	{
	}

	Fimage(const Ivec2 & size):
		Fimage(size[0], size[1])
	{
	}

	inline Type & at(const int x, const int y)
	{
		assert(0 <= x && x < mSize[0]);
		assert(0 <= y && y < mSize[1]);
		return mData[y * mSize[0] + x];
	}

	inline const Type & at(const int x, const int y) const
	{
		assert(0 <= x && x < mSize[0]);
		assert(0 <= y && y < mSize[1]);
		return mData[y * mSize[0] + x];
	}

	inline Type & at(const Ivec2 & p)
	{
		assert(0 <= p[0] && p[0] < mSize[0]);
		assert(0 <= p[1] && p[1] < mSize[1]);
		return mData[p[1] * mSize[0] + p[0]];
	}

	inline const Type & at(const Ivec2 & p) const
	{
		assert(0 <= p[0] && p[0] < mSize[0]);
		assert(0 <= p[1] && p[1] < mSize[1]);
		return mData[p[1] * mSize[0] + p[0]];
	}

	Type evalTexel(int x, int y) const
	{
		if (mWrapS == Fimage::Repeat)
			x = Math::PositiveMod(x, mSize[0]);
		else if (mWrapS == Fimage::Clamp)
			x = Math::Clamp(x, 0, mSize[0] - 1);
		else if (mWrapS == Fimage::Mirror)
			throw std::runtime_error("unimpl");
		else
			throw std::runtime_error("unimpl");

		if (mWrapT == Fimage::Repeat)
			y = Math::PositiveMod(y, mSize[1]);
		else if (mWrapT == Fimage::Clamp)
			y = Math::Clamp(y, 0, mSize[1] - 1);
		else if (mWrapT == Fimage::Mirror)
			throw std::runtime_error("unimpl");
		else
			throw std::runtime_error("unimpl");

		return this->at(x, y);
	}

	Type evalNearest(const Vec2 & uv) const
	{
		// taken from mitsuba mipmap.h
		return evalTexel(Math::FloorToInt(uv[0] * mSize[0]), Math::FloorToInt(uv[1] * mSize[1]));
	}

	Type evalBilinear(const Vec2 & uv) const
	{
		// taken from mitsuba mipmap.h
		Vec2 scaledUv = uv * Vec2(mSize) - Vec2(0.5);

		int32_t xPos = Math::FloorToInt(scaledUv[0]), yPos = Math::FloorToInt(scaledUv[1]);
		double dx1 = scaledUv[0] - xPos, dy1 = scaledUv[1] - yPos;
		double dx2 = 1.0 - dx1, dy2 = 1.0 - dy1;

		return evalTexel(xPos, yPos) * dx2 * dy2
			+ evalTexel(xPos, yPos + 1) * dx2 * dy1
			+ evalTexel(xPos + 1, yPos) * dx1 * dy2
			+ evalTexel(xPos + 1, yPos + 1) * dx1 * dy1;
	}

	Type average() const
	{
		Type result;
		forEachPixel([&](const Ivec2, const Type & color)
		{
			result += color;
		});
		result /= static_cast<double>(getNumPixels());
		return result;
	}

	void forEachPixel(const std::function<void(const Ivec2 & pixel, Type * color)> & func)
	{
		for (int y = 0;y < mSize[1];y++)
			for (int x = 0;x < mSize[0];x++)
			{
				func(Ivec2(x, y), &mData[y * mSize[0] + x]);
			}
	}

	void forEachPixel(const std::function<void(const Ivec2 & pixel, const Type & color)> & func) const
	{
		for (int y = 0; y < mSize[1]; y++)
			for (int x = 0; x < mSize[0]; x++)
			{
				func(Ivec2(x, y), mData[y * mSize[0] + x]);
			}
	}

	void forEachPixel(const std::function<void(const Vec2 & uv, Type * color)> & func)
	{
		Vec2 invResolution = 1.0f / Vec2(mSize);
		for (int y = 0;y < mSize[1];y++)
			for (int x = 0;x < mSize[0];x++)
			{
				Vec2 uv = Vec2(x + 0.5f, y + 0.5f) * invResolution;
				func(uv, &mData[y * mSize[0] + x]);
			}
	}

	void forEachPixel(const std::function<void(const Vec2 & uv, const Type & color)> & func) const
	{
		Vec2 invResolution = 1.0f / Vec2(mSize);
		for (int y = 0;y < mSize[1];y++)
			for (int x = 0;x < mSize[0];x++)
			{
				Vec2 uv = Vec2(x + 0.5f, y + 0.5f) * invResolution;
				func(uv, mData[y * mSize[0] + x]);
			}
	}

	inline const int getNumPixels() const
	{
		return mSize[0] * mSize[1];
	}

	inline Type* row(int r)
	{
		assert(0 <= r && r < mSize[1]);
		return mData[mSize[0] * r];
	}

	inline const Type* row(int r) const
	{
		assert(0 <= r && r < mSize[1]);
		return mData[mSize[0] * r];
	}

	void setZero()
	{
		std::memset(&mData[0][0], 0, sizeof(Vec3) * mSize[0] * mSize[1]);
	}

	Fimage operate(const Fimage & image, const std::function<Type(const Type & in, const Type & out)> & func) const
	{
		assert(image.mSize == this->mSize);
		Fimage result(this->mSize);
		for (int y = 0; y < this->mSize[1]; y++)
			for (int x = 0; x < this->mSize[0]; x++)
				result.at(x, y) = func(this->at(x, y), image.at(x, y));
		return result;
	}

	Fimage operator+(const Fimage & image) const { return operate(image, [](const Type & a, const Type & b) { return a + b; }); }
	Fimage operator*(const Fimage & image) const { return operate(image, [](const Type & a, const Type & b) { return a * b; }); }

	WrapMode mWrapS = WrapMode::Repeat;
	WrapMode mWrapT = WrapMode::Repeat;
	bool numValidChannels = 0;
	std::vector<Type> mData;
	Ivec2 mSize;
};

struct FimageIo
{
	enum Channel
	{
		Red = 0,
		Green = 1,
		Blue = 2
	};

	template <typename Type>
	static void LoadPfm(Fimage<Type> * result,
						int * nChannels,
						const filesystem::path &filepath,
						const std::function<Type(const Vec4 & v)> & func)
	{
		// check open file
		std::ifstream is(filepath, std::ios::binary | std::ios::in);
		if (!is.is_open()) throw std::runtime_error("FimageIo: cannot find specified PFM file");

		// read pfm file header
		int height, width;
		std::string pfm, minusOne;
		is >> pfm >> width >> height >> minusOne;
		char c = (char)is.get(); if (c == '\r') c = (char)is.get();

		// check pfm ill-conditions
		if (pfm != "PF" && pfm != "Pf") return;
		if (width < 0 || height < 0) return;

		// prepare return data
		*result = Fimage<Type>(width, height);
		int numChannels = (pfm == "PF") ? 3 : 1;
		if (nChannels) *nChannels = numChannels;

		// read pfm file
		std::vector<float> row(numChannels * width, 0.0);
		for (int y = height; y > 0; y--)
		{
			is.read((char*)(&row[0]), sizeof(float) * numChannels * width);
			for (int x = 0; x < width; x++)
			{
				if (numChannels == 1)
				{
					result->at(x, y - 1) = func(Vec4(row[x]));
				}
				else if (numChannels == 3)
				{
					result->at(x, y - 1) = func(Vec4(row[x * 3 + 0], row[x * 3 + 1], row[x * 3 + 2], 1.0));
				}
			}
		}

		// close
		is.close();
	}

	template<typename Type>
	static void LoadHdr(Fimage<Type> * result,
						int * nChannels,
						const filesystem::path & filepath,
						const std::function<Type(const Vec4 & v)> & func)
	{
		int width, height, numChannels;
		float * data = stbi_loadf(filepath.string().c_str(), &width, &height, &numChannels, 4);
		if (!data) throw std::runtime_error("FimageIo: cannot find specified HDR file");
		
		// prepare return data
		*result = Fimage<Type>(width, height);
		if (nChannels) *nChannels = numChannels;

		// read hdr file
		result->forEachPixel([&](const Ivec2 & pixel, Type * color)
		{
			int pixelIndex = pixel[1] * width + pixel[0];
			Vec4 tmpColor;
			for (int c = 0; c < 4; c++) tmpColor[c] = data[pixelIndex * 4 + c];
			*color = func(tmpColor);
		});
		stbi_image_free(data);
	}

	template<typename Type>
	static void LoadExr(Fimage<Type> * result,
						int * nChannels,
						const filesystem::path &filepath,
						const std::function<Type(const Vec4 & v)> & func)
	{
		float * out;
		int width;
		int height;
		const char * err = nullptr;

		int ret = LoadEXR(&out, &width, &height, filepath.string().c_str(), &err);

		// prepare return data
		*result = Fimage<Type>(width, height);
		if (nChannels) *nChannels = 4;

		if (ret != TINYEXR_SUCCESS)
		{
			if (err)
			{
				std::cout << "TinyExr Error : " << err << std::endl;
				FreeEXRErrorMessage(err);
			}
			throw std::runtime_error("exr loading not success");
		}

		// read exr file
		for (Uint y = 0; y < height; y++)
			for (Uint x = 0; x < width; x++)
			{
				Vec4 tmpColor(out[(y * width + x) * 4],
							  out[(y * width + x) * 4 + 1],
							  out[(y * width + x) * 4 + 2],
							  out[(y * width + x) * 4 + 3]);
				result->mData[y * width + x] = func(tmpColor);
			}

		free(out);
	}

	template <typename Type>
	static void LoadLdr(Fimage<Type> * result,
						int * nChannels,
						const filesystem::path &filepath,
						const std::function<Type(const Vec4 & v)> & func)
	{
		int width, height, numChannels;
		stbi_uc * data = stbi_load(filepath.string().c_str(), &width, &height, &numChannels, STBI_rgb_alpha);
		if (data == nullptr) throw std::runtime_error(std::string("cannot open file ") + filepath.string());

		// prepare return data
		*result = Fimage<Type>(width, height);
		if (nChannels) *nChannels = numChannels;

		// read file
		result->forEachPixel([&](const Ivec2 & pixel, Type * color) {
			int offset = (height - pixel[1] - 1) * width + pixel[0];
			Vec4 tmpColor;
			for (Uint c = 0; c < 4; c++)
			{
				tmpColor[c] = (c < numChannels) ? pow(static_cast<double>(data[offset * 4 + c]) / 255.0, 2.2) : 1.0;
			}
			*color = func(tmpColor);
		});

		stbi_image_free(data);
	}

	template <typename Type>
	static void LoadAny(Fimage<Type> * result,
						int * nChannels,
						const filesystem::path & filepath,
						const std::function<Type(const Vec4 & v)> & func)
	{
		if (!filepath.has_extension()) throw std::runtime_error("unknown file type");
		filesystem::path extension = filepath.extension();
		if (extension == ".pfm")
			LoadPfm<Type>(result, nChannels, filepath, func);
		else if (extension == ".hdr")
			LoadHdr<Type>(result, nChannels, filepath, func);
		else if (extension == ".exr")
			LoadExr<Type>(result, nChannels, filepath, func);
		else if (extension == ".png" || extension == ".jpg" || extension == ".bmp" || extension == ".tga")
			LoadLdr<Type>(result, nChannels, filepath, func);
		else
			throw std::runtime_error("unsupported file format");
	}

	static Fimage<Spectrum> Load(int * nChannels, const filesystem::path & filepath)
	{
		auto v = [&](const Vec4 & vec) { return SpectrumConvertUtil::SpectrumFromRgb(vec[0], vec[1], vec[2]); };
		Fimage<Spectrum> result;
		LoadAny<Spectrum>(&result, nChannels, filepath, v);
		return result;
	}

	static shared_ptr<Fimage<Spectrum>> LoadPtr(int * nChannels, const filesystem::path & filepath)
	{
		auto v = [&](const Vec4 & vec) { return SpectrumConvertUtil::SpectrumFromRgb(vec[0], vec[1], vec[2]); };
		shared_ptr<Fimage<Spectrum>> result = make_shared<Fimage<Spectrum>>();
		LoadAny<Spectrum>(result.get(), nChannels, filepath, v);
		return result;
	}

	static Fimage<double> Load(int * nChannels, const filesystem::path & filepath, const int channel)
	{
		auto v = [&](const Vec4 & vec) { return vec[channel]; };
		Fimage<double> result;
		LoadAny<double>(&result, nChannels, filepath, v);
		return result;
	}

	static shared_ptr<Fimage<double>> LoadPtr(int * nChannels, const filesystem::path & filepath, const int channel)
	{
		auto v = [&](const Vec4 & vec) { return vec[channel]; };
		shared_ptr<Fimage<double>> result = make_shared<Fimage<double>>();
		LoadAny<double>(result.get(), nChannels, filepath, v);
		return result;
	}

	static Fimage<double> Load(const filesystem::path & filepath, const int channel) { return Load(nullptr, filepath, channel); }
	static Fimage<Spectrum> Load(const filesystem::path & filepath) { return Load(nullptr, filepath); }
	static shared_ptr<Fimage<double>> LoadPtr(const filesystem::path & filepath, const int channel) { return LoadPtr(nullptr, filepath, channel); }
	static shared_ptr<Fimage<Spectrum>> LoadPtr(const filesystem::path & filepath) { return LoadPtr(nullptr, filepath); }

	static void SavePfmColor(const Fimage<Spectrum> &fimage,
							 const filesystem::path &filepath,
							 const std::function<Vec4(const Spectrum & v)> & func)
	{
		int width = fimage.mSize[0];
		int height = fimage.mSize[1];

		// write PFM header
		std::ofstream os(filepath, std::ios::binary);
		if (!os.is_open()) throw std::runtime_error("cannot open file");
		os << "PF" << std::endl;
		os << width << " " << height << std::endl;
		os << "-1" << std::endl;

		// write data flip row
		for (int y = height;y >= 1;y--)
			for (int x = 0; x < width; x++)
			{
				const Vec4 v = func(fimage.at(x, y - 1));
				for (int c = 0; c < 3; c++)
				{
					float vf = float(v[c]);
					os.write((char*)(&vf), sizeof(float));
				}
			}

		os.close();
	}

	static void SavePfmMono(const Fimage<double> &fimage, const filesystem::path &filepath)
	{
		int width = fimage.mSize[0];
		int height = fimage.mSize[1];

		// write PFM header
		std::ofstream os(filepath, std::ios::binary);
		if (!os.is_open()) throw std::runtime_error("cannot open file");
		os << "Pf" << std::endl;
		os << width << " " << height << std::endl;
		os << "-1" << std::endl;

		// write data flip row
		for (int y = height; y >= 1; y--)
			for (int x = 0; x < width; x++)
			{
				float vf = float(fimage.at(x, y - 1));
				os.write((char*)(&vf), sizeof(float));
			}

		os.close();
	}

	template <typename Type>
	static void SaveHdr(const Fimage<Type> & fimage,
						const filesystem::path &filepath,
						const std::function<Vec4(const Type & v)> & func)
	{
		float * data = new float[fimage.mSize[0] * fimage.mSize[1] * 3];
		for (int row = 0;row < fimage.mSize[1];row++)
			for (int col = 0; col < fimage.mSize[0]; col++)
			{
				const Vec4 v = func(fimage.at(col, row));
				for (int i = 0; i < 3; i++)
				{
					data[(col + row * fimage.mSize[0]) * 3 + i] = static_cast<float>(v[i]);
				}
			}

		stbi_write_hdr(filepath.string().c_str(),
					   static_cast<int>(fimage.mSize[0]),
					   static_cast<int>(fimage.mSize[1]),
					   3,
					   data);
		delete [] data;
	}

	template <typename Type>
	static void SavePng(const Fimage<Type> & fimage,
						const filesystem::path &filepath,
						const std::function<Vec4(const Type & v)> & func)
	{
		char * data = new char[fimage.mSize[0] * fimage.mSize[1] * 3];
		for (int row = 0;row < fimage.mSize[1];row++)
			for (int col = 0; col < fimage.mSize[0]; col++)
			{
				const Vec4 v = Math::Clamp(func(fimage.at(col, row)), Vec4(0.0), Vec4(1.0));
				for (Uint i = 0; i < 3; i++)
				{
					double p = std::pow(v[i], double(1 / 2.2));
					p = std::min(p * 255.9999, 255.0);
					data[(col + row * fimage.mSize[0]) * 3 + i] = static_cast<char>(p);
				}
			}

		stbi_write_png(filepath.string().c_str(),
			static_cast<int>(fimage.mSize[0]),
			static_cast<int>(fimage.mSize[1]),
			3,
			data,
			static_cast<int>(fimage.mSize[0] * 3));
		delete [] data;
	}

	template <typename Type>
	static void SaveAny(const Fimage<Type> & fimage,
						const filesystem::path & filepath,
						const std::function<Vec4(const Type & v)> & func)
	{
		if (!filepath.has_extension()) throw std::runtime_error("unknown file type");
		filesystem::path extension = filepath.extension();

		if (extension == ".hdr")
			SaveHdr<Type>(fimage, filepath, func);
		else if (extension == ".png")
			SavePng<Type>(fimage, filepath, func);
		else
			throw std::runtime_error("unsupported file format");
	}

	static void Save(const Fimage<double> & fimage, const filesystem::path & filepath)
	{
		// pfm for mono is special
		filesystem::path extension = filepath.extension();
		if (extension == ".pfm")
		{
			SavePfmMono(fimage, filepath);
			return;
		}

		// otherwise just follow default procedure
		auto func = [&](const double & v) { return Vec4(v, v, v, 1.0); };
		SaveAny<double>(fimage, filepath, func);
	}

	static void Save(const Fimage<RgbSpectrum> & fimage, const filesystem::path & filepath)
	{
		auto func = [&](const RgbSpectrum & v) { Vec3 t = SpectrumConvertUtil::Vec3FromSpectrum(v); return Vec4(t[0], t[1], t[2], 1.0); };

		// pfm for mono is special
		filesystem::path extension = filepath.extension();
		if (extension == ".pfm")
		{
			SavePfmColor(fimage, filepath, func);
			return;
		}

		// otherwise just follow default procedure
		SaveAny<RgbSpectrum>(fimage, filepath, func);
	}
};