#pragma once

#include "common/wurst.h"
#include "common/cdftable.h"
#include "accel/embreeaccel.h"
#include "light/envmap.h"

struct Scene
{
	Scene(const shared_ptr<const EmbreeAccel> & accel,
		  const std::vector<shared_ptr<const Geometry>> & geometries,
		  const std::vector<shared_ptr<const Geometry>> & lightGeometries,
		  shared_ptr<EnvmapLight> envmap,
		  const double planeZ = 0.0):
		mAccel(accel),
		mGeometries(geometries),
		mLightGeometries(lightGeometries),
		mEnvmapLight(envmap),
		mPlaneZ(planeZ)
	{
		if (envmap)
		{
			Vec3 center; double radius = 0.0;
			accel->boundingSphere(&center, &radius);
			envmap->updateSceneRadius(center, radius);
		}
	}

	Scene(const shared_ptr<const EmbreeAccel> & accel,
		  const std::vector<shared_ptr<const Geometry>> & geometries,
		  const std::vector<shared_ptr<const Geometry>> & lightGeometries,
		  const double planeZ = 0.0):
		Scene(accel, geometries, lightGeometries, nullptr, planeZ)
	{
	}

	Scene(const shared_ptr<const EmbreeAccel> & accel,
		  const std::vector<shared_ptr<const Geometry>> & geometries,
		  shared_ptr<EnvmapLight> envmap,
		  const double planeZ = 0.0):
		mAccel(accel),
		mGeometries(geometries),
		mEnvmapLight(envmap),
		mPlaneZ(planeZ)
	{
		if (envmap)
		{
			Vec3 center; double radius = 0.0;
			accel->boundingSphere(&center, &radius);
			envmap->updateSceneRadius(center, radius);
		}
	}

	Scene(const shared_ptr<const EmbreeAccel> & accel,
		  const std::vector<shared_ptr<const Geometry>> & geometries,
		  const double planeZ = 0.0):
		Scene(accel, geometries, nullptr, planeZ)
	{
	}


	Bound2 bound2() const
	{
		Bound2 result;
		for (const shared_ptr<const Geometry> & geometry : mGeometries)
		{
			result = Bound2::Union(result, geometry->bound2());
		}
		return result;
	}

	Bound3 bound() const
	{
		return mAccel->bound();
	}

	void boundingSphere(Vec3 * position, double * radius) const
	{
		mAccel->boundingSphere(position, radius);
	}

	shared_ptr<Vertex> intersect(Ray3 * ray) const
	{
		shared_ptr<Vertex> v = mAccel->intersect(ray);
		if (v) { return v; }
		if (mEnvmapLight) return mEnvmapLight->intersect(ray);
		return nullptr;
	}

	bool isIntersect(const Ray3 & ray) const
	{
		return mAccel->isIntersect(ray);
	}

	Spectrum transmittance(const Vertex & p1, const Vertex & p2) const
	{
		assert(!(p1.hasType(VertexType::Envmap) && p2.hasType(VertexType::Envmap)));
		if (p1.hasType(VertexType::Envmap))
		{
			const double vis = isIntersect(Ray3(p2.mPosition, p1.mPosition)) ? 0.0 : 1.0;
			return Spectrum(vis);
		}
		else if (p2.hasType(VertexType::Envmap))
		{
			const double vis = isIntersect(Ray3(p1.mPosition, p2.mPosition)) ? 0.0 : 1.0;
			return Spectrum(vis);
		}
		else
		{
			const Vec3 uv = p2.mPosition - p1.mPosition;
			const double length = Math::Length(uv);

			// handle volumetric transmission
			Spectrum tr(1.0);
			const Medium * uMedium = p1.getMedium(uv);
			const Medium * vMedium = p2.getMedium(-uv);

#if 0
			if (uMedium == vMedium)
			{
				const Medium * medium = uMedium;
				if (medium) tr = medium->evalTransmittance(length);
			}
			else
			{
				tr = Spectrum(0.0);
			}
#else
			const Medium * medium = nullptr;
			if (vMedium != nullptr);
			{
				medium = vMedium;
			}
			if (uMedium != nullptr)
			{
				medium = uMedium;
			}
			if (medium) tr = medium->evalTransmittance(length);
#endif

			// handle visibility
			const Vec3 nuv = uv / length;
			const double vis = isIntersect(Ray3(p1.mPosition, nuv, Ray3::Epsilon, length - Ray3::Epsilon)) ? 0.0 : 1.0;

			return tr * vis;
		}
	}

	shared_ptr<const EmbreeAccel> mAccel = nullptr;
	shared_ptr<const EnvmapLight> mEnvmapLight = nullptr;
	std::vector<shared_ptr<const Geometry>> mGeometries;
	std::vector<shared_ptr<const Geometry>> mLightGeometries;
	double mPlaneZ = 0.0; // required for 2d global illumination
};
