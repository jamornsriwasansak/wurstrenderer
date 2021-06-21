#pragma once

#include "common/wurst.h"

#include "common/bsdf.h"
#include "common/camera.h"
#include "common/light.h"
#include "common/medium.h"
#include "common/phase.h"
#include "common/subpathsampler.h"
#include "common/util/file.h"
#include "common/sphericalrepresentation.h"

// vertex type start from light source and end at camera
enum class VertexType : uint8_t
{
	LightImaginary = 1 << 0,
	Light = 1 << 1,
	Envmap = 1 << 2,
	Surface = 1 << 3,
	Medium = 1 << 4,
	Camera = 1 << 5,
	CameraImaginary = 1 << 6
};

DECLARE_ENUM_OPERATORS(VertexType);

struct Vertex
{
	static void AssertOrder(const Vertex & u, const Vertex & v)
	{
		assert(!u.hasType(VertexType::CameraImaginary));
		if (u.hasType(VertexType::LightImaginary)) { assert(v.hasType(VertexType::Light | VertexType::Envmap)); }
		if (u.hasType(VertexType::Light | VertexType::Envmap)) { assert(v.hasType(VertexType::Surface | VertexType::Medium | VertexType::Camera)); }
		if (u.hasType(VertexType::Surface)) { assert(v.hasType(VertexType::Surface | VertexType::Medium | VertexType::Camera)); }
		if (u.hasType(VertexType::Medium)) { assert(v.hasType(VertexType::Surface | VertexType::Medium | VertexType::Camera)); }
		if (u.hasType(VertexType::Camera)) { assert(v.hasType(VertexType::CameraImaginary)); }
	}

	static void AssertOrder(const Vertex & u, const Vertex & mid, const Vertex & v)
	{
		const Vertex & p1 = (u.getVertexType() <= v.getVertexType()) ? u : v;
		const Vertex & p2 = (u.getVertexType() <= v.getVertexType()) ? v : u;
		AssertOrder(p1, mid); AssertOrder(mid, p2);
	}

	static double ConvertToAreaFromSolidAngle(const Vertex & from, const Vertex & to)
	{
		assert(!from.hasType(VertexType::CameraImaginary | VertexType::LightImaginary));
		assert(!to.hasType(VertexType::CameraImaginary | VertexType::LightImaginary));
        assert(!(from.hasType(VertexType::Envmap) && to.hasType(VertexType::Envmap)));

		double result = 1.0;
		Vec3 uv = to - from;
		double length2 = Math::Length2(uv);
		Vec3 nuv = uv / std::sqrt(length2);
		if (!to.hasType(VertexType::Medium | VertexType::Envmap)) { result *= Math::AbsDot(to.mGeometryNormal, nuv); }
		result /= length2;
		return result;
	}

	static double GeometryTerm(const Vertex & u, const Vertex & v)
	{
		assert(!u.hasType(VertexType::CameraImaginary | VertexType::LightImaginary));
		assert(!v.hasType(VertexType::CameraImaginary | VertexType::LightImaginary));
		assert(!(u.hasType(VertexType::Envmap) && v.hasType(VertexType::Envmap)));
		double result = 1.0;
		Vec3 uv = v - u;
		double length2 = Math::Length2(uv);
		Vec3 nuv = uv / std::sqrt(length2);
		// compute cos cos / d^2
		if (!u.hasType(VertexType::Medium | VertexType::Envmap)) { result *= Math::AbsDot(u.mGeometryNormal, nuv); }
		if (!v.hasType(VertexType::Medium | VertexType::Envmap)) { result *= Math::AbsDot(v.mGeometryNormal, nuv); }
		result /= length2;
		return result;
	}

	template<typename V>
	static shared_ptr<V> SafeCast(shared_ptr<Vertex> v)
	{
		assert(v->getVertexType() == V().getVertexType());
		return static_pointer_cast<V>(v);
	}

	Vertex(): mPosition(Vec3(0.0)) {}

	~Vertex() {}

	bool hasType(VertexType vt) const { return (this->getVertexType() & vt) != 0; } 

	const Medium * getMedium(const Vec3 & direction) const { return Math::IsFrontface(direction, mGeometryNormal) ? mMediumFrontface : mMediumBackface; }

	Spectrum scatter(const Vertex & v1, const Vertex & v2) const { return scatterImpl(nullptr, v1, v2); }

	Spectrum scatter(const Vec3 & wi, const Vec3 & wo) const { return scatterImpl(nullptr, wi, wo); }

	Spectrum scatter(Vec2 * raster, const Vertex & v1, const Vertex & v2) const { return scatterImpl(raster, v1, v2); }

	Spectrum scatter(Vec2 * raster, const Vec3 & wi, const Vec3 & wo) const { throw scatterImpl(raster, wi, wo); }

	virtual double pdfA(const Vertex & v1, const Vertex & v2) const { return pdfW(v1, v2) * Vertex::ConvertToAreaFromSolidAngle(*this, v2); };

	// common

	virtual VertexType getVertexType() const = 0;

	// 3d case, pdf

	virtual double pdfW(const Vertex & v1, const Vertex & v2) const { throw std::runtime_error("unimpl"); return 0.0; }

	virtual double pdfW(const Vec3 & wi, const Vec3 & wo) const { throw std::runtime_error("unimpl"); return 0.0; }

	// 3d case, scatter

	virtual Spectrum scatterImpl(Vec2 * raster, const Vertex & v1, const Vertex & v2) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	virtual Spectrum scatterImpl(Vec2 * raster, const Vec3 & wi, const Vec3 & wo) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	virtual Spectrum source() const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

    Vec3 operator-(const Vertex & from) const
    {
	    const Vertex & to = *this;
        assert(!from.hasType(VertexType::CameraImaginary | VertexType::LightImaginary));
        assert(!to.hasType(VertexType::CameraImaginary | VertexType::LightImaginary));
        assert(!(from.hasType(VertexType::Envmap) && to.hasType(VertexType::Envmap)));

		// the position of envmap is basically infinitely far from the scene
		// therefore dir_to_env = Normalize(x_env - x_any) = x_env which match what we store in the envmap mPosition exactly
        if (from.hasType(VertexType::Envmap)) { return -from.mPosition; }
        if (to.hasType(VertexType::Envmap)) { return to.mPosition; }

        return to.mPosition - from.mPosition;
    }

	// 2d case, pdf2

	// 2d case, scatter2

	Spectrum scatter2(const Vertex & v1, const Vertex & v2) const { return scatter2(nullptr, v1, v2); }

	virtual Spectrum scatter2(Vec2 * raster, const Vertex & v1, const Vertex & v2) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	virtual Spectrum scatter2(Vec2 * raster, const Vec2 & v1, const Vec2 & v2) const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }

	virtual Spectrum source2() const { throw std::runtime_error("unimpl"); return Spectrum(0.0); }


	// contribution (f(X)/p(X)) of the path until this vertex
	Spectrum mPathContrib = Spectrum(0.0);

	const SubpathSampler * mSubpathSampler;

	// I have committed a crime.

	SUnion
	{
		Vec3 mPosition;
		Vec2 mPosition2;
	};

	SUnion
	{
		Vec3 mGeometryNormal;
		Vec2 mGeometryNormal2;
	};
	
	SUnion
	{
		CoordFrame3 mCoordFrame;
		CoordFrame2 mCoordFrame2;
	};

	const Medium * mMediumBackface = nullptr;
	const Medium * mMediumFrontface = nullptr;
};

struct LightImaginaryVertex : public Vertex
{
	VertexType getVertexType() const override { return VertexType::LightImaginary; }

	LightImaginaryVertex()
	{
	}
};

struct LightVertex : public Vertex
{
	VertexType getVertexType() const override { return VertexType::Light; }

	LightVertex()
	{
	}

	// 3d case

	double pdfW(const Vertex & v1, const Vertex & v2) const override
	{
		return pdfW(Vec3(0.0), Math::Normalize(v2 - *this));
	}

	double pdfW(const Vec3 & wi, const Vec3 & wo) const override
	{
		return mLight->evalPdfLe1(*this, mCoordFrame.toLocal(wo));
	}

	using Vertex::scatter;

	Spectrum scatterImpl(Vec2 *, const Vertex & v1, const Vertex & v2) const override
	{
		const Vertex & p1 = (v1.getVertexType() <= v2.getVertexType()) ? v1 : v2;
		const Vertex & p2 = (v1.getVertexType() <= v2.getVertexType()) ? v2 : v1;
		AssertOrder(p1, *this); AssertOrder(*this, p2);
		return scatterImpl(nullptr, Vec3(0.0), Math::Normalize(p2 - *this));
	}

	Spectrum scatterImpl(Vec2 * , const Vec3 & , const Vec3 & wo) const override
	{
		return mLight->evalLe1(*this, mCoordFrame.toLocal(wo));
	}

	Spectrum source() const override
	{
		return mLight->evalLe0(*this);
	}

	// 2d case

	Spectrum scatter2(Vec2 *, const Vertex & v1, const Vertex & v2) const override
	{
		const Vertex & p1 = (v1.getVertexType() <= v2.getVertexType()) ? v1 : v2;
		const Vertex & p2 = (v1.getVertexType() <= v2.getVertexType()) ? v2 : v1;
		AssertOrder(p1, *this); AssertOrder(*this, p2);
		return scatter2(nullptr, Vec2(0.0), Math::Normalize(p2.mPosition2 - mPosition2));
	}

	Spectrum scatter2(Vec2 * , const Vec2 & , const Vec2 & wo) const override
	{
		return mLight->evalLe12(*this, mCoordFrame2.toLocal(wo));
	}

	Spectrum source2() const override
	{
		return mLight->evalLe02(*this);
	}

	// members

	Vec2 mTextureCoord;
	const Light * mLight = nullptr;
	const Geometry * mGeometry = nullptr;
};

struct EnvmapVertex : public LightVertex
{
	VertexType getVertexType() const override { return VertexType::Envmap | VertexType::Light; }

	EnvmapVertex()
	{
	}

	// 3d case

	double pdfW(const Vertex & v1, const Vertex & v2) const override { throw std::runtime_error("shouldn't be called"); return 0.0; }

	double pdfW(const Vec3 & wi, const Vec3 & wo) const override { throw std::runtime_error("shouldn't be called"); return 0.0; }

	double pdfA(const Vertex & v1, const Vertex & v2) const override
	{
		assert(!v2.hasType(VertexType::Envmap));
		// first compute the probably of sampling the disk
		double diskPdfA = mLight->evalPdfLe1(*this, Vec3(0.0));
		return diskPdfA * Vertex::ConvertToAreaFromSolidAngle(*this, v2);
	}

	using Vertex::scatter;

	Spectrum scatterImpl(Vec2 *, const Vertex & v1, const Vertex & v2) const override
	{
		const Vertex & p1 = (v1.getVertexType() <= v2.getVertexType()) ? v1 : v2;
		const Vertex & p2 = (v1.getVertexType() <= v2.getVertexType()) ? v2 : v1;
		AssertOrder(p1, *this); AssertOrder(*this, p2);
		return scatterImpl(nullptr, Vec3(0.0), Vec3(0.0));
	}

	Spectrum scatterImpl(Vec2 * , const Vec3 & , const Vec3 & wo) const override
	{
		return mLight->evalLe1(*this, mPosition);
	}

	Spectrum source() const override
	{
		return mLight->evalLe0(*this);
	}
};

struct SurfaceVertex : public Vertex
{
	VertexType getVertexType() const override { return VertexType::Surface; }

	SurfaceVertex()
	{
	}

	// 3d case

	double pdfW(const Vertex & v1, const Vertex & v2) const override
	{
		return pdfW(Math::Normalize(v1 - *this), Math::Normalize(v2 - *this));
	}

	double pdfW(const Vec3 & wi, const Vec3 & wo) const override
	{
		return mBsdf->evalPdfW(*this, mCoordFrame.toLocal(wi), mCoordFrame.toLocal(wo));
	}

	using Vertex::scatter;

	Spectrum scatterImpl(Vec2 *, const Vertex & v1, const Vertex & v2) const override
	{
		AssertOrder(v1, *this, v2);
		return scatterImpl(nullptr, Math::Normalize(v1 - *this), Math::Normalize(v2 - *this));
	}

	Spectrum scatterImpl(Vec2 * , const Vec3 & wi, const Vec3 & wo) const override
	{
		return mBsdf->evalBsdf(*this, mCoordFrame.toLocal(wi), mCoordFrame.toLocal(wo));
	}

	// 2d case

	// TODO::

	// members

	Vec2 mTextureCoord;
	const Bsdf * mBsdf = nullptr;
	const Geometry * mGeometry = nullptr;
};

struct MediumVertex : public Vertex
{
	VertexType getVertexType() const override { return VertexType::Medium; }

	double pdfW(const Vertex & v1, const Vertex & v2) const override
	{
		return pdfW(Math::Normalize(v1 - *this), Math::Normalize(v2 - *this));
	}

	double pdfW(const Vec3 & wi, const Vec3 & wo) const override
	{
		return mPhaseFunction->evalPdfW(*this, mCoordFrame.toLocal(wi), mCoordFrame.toLocal(wo));
	}

	using Vertex::scatter;

	Spectrum scatterImpl(Vec2 *, const Vertex & v1, const Vertex & v2) const override
	{
		AssertOrder(v1, *this, v2);
		return scatterImpl(nullptr, Math::Normalize(v1 - *this), Math::Normalize(v2 - *this));
	}

	Spectrum scatterImpl(Vec2 * , const Vec3 & wi, const Vec3 & wo) const override
	{
		return mPhaseFunction->evalPhase(*this, mCoordFrame.toLocal(wi), mCoordFrame.toLocal(wo));
	}

	// 2d case

	// TODO::

	// members

	const PhaseFunction * mPhaseFunction;
};

struct CameraVertex : public Vertex
{
	VertexType getVertexType() const override { return VertexType::Camera; }

	double pdfW(const Vertex & v1, const Vertex & v2) const override
	{
		assert(v1.hasType(VertexType::CameraImaginary));
		return mSubpathSampler->evalCameraPdfW(*this, Math::Normalize(v2 - *this));
	}

	double pdfW(const Vec3 & , const Vec3 & wo) const override
	{
		return mSubpathSampler->evalCameraPdfW(*this, wo);
	}

	using Vertex::scatter;

	Spectrum scatterImpl(Vec2 * raster, const Vertex & v1, const Vertex & v2) const override
	{
		const Vertex & p1 = (v1.getVertexType() <= v2.getVertexType()) ? v1 : v2;
		const Vertex & p2 = (v1.getVertexType() <= v2.getVertexType()) ? v2 : v1;
		AssertOrder(p1, *this); AssertOrder(*this, p2);
		return scatterImpl(raster, Math::Normalize(p1 - *this), Vec3(0.0));
	}

	Spectrum scatterImpl(Vec2 * raster, const Vec3 & wi, const Vec3 &) const override
	{
		return mCamera->evalWe1(raster, *this, wi);
	}

	Spectrum source() const override
	{
		return mCamera->evalWe0(*this);
	}

	SUnion
	{
		const Camera * mCamera;
		const Camera2 * mCamera2;
	};
};

struct CameraImaginaryVertex : public Vertex
{
	using Vertex::scatter;

	VertexType getVertexType() const override { return VertexType::CameraImaginary; }
};
