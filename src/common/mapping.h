#pragma once

#include <iostream>
#include <array>
#include "math/math.h"

struct FiboInterpolateInfo
{
	int index;
	double weight;
};

struct Mapping
{
	static inline Vec2 ArcFromSegment(const double sample, const double halfAngleMax)
	{
		return World2FromAngle(2.0 * halfAngleMax * sample - halfAngleMax);
	}

	// taken from https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
	static inline Vec2 BarycentricFromSquare(const Vec2 & sample)
	{
		const double sqrtSample0 = std::sqrt(sample[0]);
		double b1 = (sqrtSample0 * (1.0 - sample[1]));
		double b2 = (sample[1] * sqrtSample0);
		return Vec2(b1, b2);
	}

	static inline Vec2 CircleFromSegment(const double sample)
	{
		const double theta = 2.0 * Math::Pi * sample;
		return World2FromAngle(theta);
	}

	static inline Vec2 CosineWeightedHemicircleFromSegment(const double & sample)
	{
		const double cosTheta = 1.0 - 2.0 * sample;
		const double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
		return Vec2(cosTheta, sinTheta);
	}

	static inline Vec3 CosineWeightedHemisphereFromSquare(const Vec2 & sample)
	{
		const double sinPhi = std::sqrt(1.0 - sample[0]);
		const double theta = Math::Pi * 2.0 * sample[1]; // theta
		const double cosTheta = std::cos(theta);
		const double sinTheta = std::sin(theta);
		return Vec3(cosTheta * sinPhi, std::sqrt(sample[0]), sinTheta * sinPhi);
	}

	// taken from Dave Cline from Peter Shirley from http://psgraphics.blogspot.jp/2011/01/improved-code-for-concentric-map.html
	// and https://github.com/mmp/pbrt-v3/blob/7095acfd8331cdefb03fe0bcae32c2fc9bd95980/src/core/sampling.cpp
	static inline Vec2 DiskFromSquare(const Vec2 & sample)
	{
		Vec2 ab = sample * (double)2.0 - Vec2(1.0);

		if (ab[0] == 0 && ab[1] == 0) { return Vec2(0.0); }

		Vec2 ab2 = ab * ab;
		double phi, r;
		if (ab2[0] > ab2[1]) { // use squares instead of absolute values
			r = ab[0];
			phi = (Math::Pi / 4.0) * (ab[1] / ab[0]);
		}
		else {
			r = ab[1];
			phi = (Math::Pi / 2.0) - (Math::Pi / 4.0) * (ab[0] / ab[1]);
		}
		double cosphi = std::cos(phi);
		double sinphi = std::sin(phi);

		return r * Vec2(cosphi, sinphi);
	}

	static inline Vec2 HemicircleFromSegment(const double sample)
	{
		const double theta = Math::Pi * (sample - 0.5);
		return World2FromAngle(theta);
	}

	static inline Vec3 HemisphereFromSquare(const Vec2 & sample)
	{
		// z = [0, 1]
		const double sinPhi = std::sqrt(1.0 - sample[0] * sample[0]);
		const double theta = 2.0 * Math::Pi * sample[1]; // theta = [0, 2pi)
		const double cosTheta = std::cos(theta);
		const double sinTheta = std::sin(theta);
		return Vec3(cosTheta * sinPhi, sample[0], sinTheta * sinPhi);
	}

	static inline Vec2 PanoramaFromWorld(const Vec3 & dir)
	{
		double u = std::atan2(-dir[2], -dir[0]) * Math::InvPi * 0.5f + 0.5f;
		double v = std::acos(dir[1]) * Math::InvPi;
		return Vec2(u, v);
	}

	// James Arvo 1992, Fast random rotation matrices
	// not truly uniform. do not use unless specified.
	static inline CoordFrame3 RotationMatrixFromCubeArvo(const Vec3 & sample)
	{
		const double theta = 2.0 * Math::Pi * sample[0];
		const double cosTheta = std::cos(theta);
		const double sinTheta = std::sin(theta);
		const double phi = 2.0 * Math::Pi * sample[1];
		const double cosPhi = std::cos(phi);
		const double sinPhi = std::sin(phi);
		const double z = sample[2];
		const double sqrtZ = std::sqrt(z);
		const double sqrtOneMinusZ = std::sqrt(1.0 - z);

		const double v00 = 2.0 * cosPhi * cosPhi * z;
		const double v01 = 2.0 * sinPhi * cosPhi * z;
		const double v02 = 2.0 * cosPhi * std::sqrt(z - z * z);
		const double v11 = 2.0 * sinPhi * sinPhi * z;
		const double v12 = 2.0 * sinPhi * std::sqrt(z - z * z);
		const double v22 = 2.0 - 2.0 * z;

		// rotation about z axis
		Matrix4 aaa(v00 - 1,	v01,		v02,		0.0, 
					v01,		v11 - 1,	v12,		0.0,
					v02,		v12,		v22 - 1,	0.0,
					0.0, 0.0, 0.0, 1.0);

		// rotation z axis to random orientation
		Matrix4 bbb(cosTheta,	-sinTheta,	0.0,	0.0,
					sinTheta,	cosTheta,	0.0,	0.0,
					0.0, 0.0, 1.0, 0.0,
					0.0, 0.0, 0.0, 1.0);

		Matrix4 m = aaa * bbb;
		CoordFrame3 result;
		result.x = Vec3(m.cols[0]);
		result.y = Vec3(m.cols[1]);
		result.z = Vec3(m.cols[2]);
		return result;
	}

	// Ken Shoemaker 1992, Uniform random rotation matrix
	static inline CoordFrame3 RotationMatrixFromCube(const Vec3 & sample)
	{
		const double theta = 2.0 * Math::Pi * sample[0];
		const double cosTheta = std::cos(theta);
		const double sinTheta = std::sin(theta);
		const double phi = 2.0 * Math::Pi * sample[1];
		const double cosPhi = std::cos(phi);
		const double sinPhi = std::sin(phi);
		const double z = sample[2];
		const double sqrtZ = std::sqrt(z);
		const double sqrtOneMinusZ = std::sqrt(1.0 - z);

		Matrix4 m = Matrix4::FromQuaternion(sqrtZ * cosPhi, sqrtOneMinusZ * sinTheta, sqrtOneMinusZ * cosTheta, sqrtZ * sinPhi);
		CoordFrame3 result;
		result.x = Vec3(m.cols[0]);
		result.y = Vec3(m.cols[1]);
		result.z = Vec3(m.cols[2]);
		return result;
	}

	// Total Compendium pg. 19 (34)
	static inline Vec3 SolidAngleFromSquare(const Vec2 & sample, const double halfAngleMax)
	{
		const double phi = 2.0 * Math::Pi * sample[0];
		const double y = 1.0 - sample[1] * (1.0 - std::cos(halfAngleMax));
		const double l = std::sqrt(1.0 - y * y);
		const double cosphi = std::cos(phi);
		const double sinphi = std::sin(phi);
		return Vec3(cosphi * l, y, sinphi * l);
	}

	static inline Vec2 SphericalFromWorld(const Vec3 & pos)
	{
		const double theta = std::atan2(pos[2], pos[0]);
		const double numerator = std::sqrt(pos[0] * pos[0] + pos[2] * pos[2]);
		const double phi = std::atan2(numerator, pos[1]);
		return Vec2(theta, phi);
	}

	static inline Vec3 SphereFromSquare(const Vec2 & sample)
	{
		const double z = 1.0 - 2.0 * sample[1];
		const double r = std::sqrt(sample[1] * (1.0 - sample[1]));
		const double phi = 2.0 * Math::Pi * sample[0]; // theta = [0, 2pi)
		const double cosphi = std::cos(phi);
		const double sinphi = std::sin(phi);
		return Vec3(2.0 * cosphi * r, 2.0 * sinphi * r, z);
	}

	// Spherical Fibonacci Mapping by Benjamin Keinert et al. section 4.
	static inline int SphFiboIndexFromWorld(const Vec3 & p, const int nSamples)
	{
		assert(nSamples > 0);
		const auto madfrac = [](const double a, const double b) { return a * b - std::floor(a * b); };

		const double invNSamples = 1.0 / double(nSamples);
		const double GoldenPHI = (std::sqrt(5.0) + 1.0) / 2.0;
		const double theta = std::min(std::atan2(p[2], p[0]), Math::Pi);
		const double cosPhi = p[1];
		const double k = std::max(2.0, std::floor(std::log(nSamples * Math::Pi * std::sqrt(5.0) * (1.0 - cosPhi * cosPhi)) / std::log(GoldenPHI * GoldenPHI)));
		const double Fk = std::pow(GoldenPHI, k) / std::sqrt(5.0);
		const double F0 = std::round(Fk);
		const double F1 = std::round(Fk * GoldenPHI);
		const Vec2 BRow0(2.0 * Math::Pi * madfrac(F0 + 1.0, GoldenPHI - 1.0) - 2.0 * Math::Pi * (GoldenPHI - 1.0),
					  2.0 * Math::Pi * madfrac(F1 + 1.0, GoldenPHI - 1.0) - 2.0 * Math::Pi * (GoldenPHI - 1.0));
		const Vec2 BRow1(-2.0 * F0 * invNSamples,
						 -2.0 * F1 * invNSamples);
		const double invDetB = 1.0 / (BRow0[0] * BRow1[1] - BRow0[1] * BRow1[0]);
		const Vec2 invBRow0(BRow1[1] * invDetB, -BRow0[1] * invDetB);
		const Vec2 invBRow1(-BRow1[0] * invDetB, BRow0[0] * invDetB);
		const Vec2 cc(theta, cosPhi - (1.0 - 1.0 * invNSamples));
		const Vec2 c = Math::Floor(Vec2(Math::Dot(cc, invBRow0), Math::Dot(cc, invBRow1)));
		double d = std::numeric_limits<double>::max();
		double j = 0;
		for (int s = 0;s < 4;s++)
		{
			const Vec2 uv(double(s % 2), double(s / 2));
			double qCosPhi = Math::Dot(BRow1, uv + c) + (1 - invNSamples);
			qCosPhi = Math::Clamp(qCosPhi, -1.0, +1.0) * 2.0 - qCosPhi;

			const double i = std::floor(nSamples * 0.5 - qCosPhi * nSamples * 0.5);
			const double qTheta = 2 * Math::Pi * madfrac(i, GoldenPHI - 1);
			qCosPhi = 1 - (2 * i + 1) * invNSamples;
			const double qSinPhi = std::sqrt(1 - qCosPhi * qCosPhi);

			const Vec3 q(std::cos(qTheta) * qSinPhi, qCosPhi, std::sin(qTheta) * qSinPhi);
			const double squaredDistance = Math::Distance2(q, p);
			if (squaredDistance < d)
			{
				d = squaredDistance;
				j = i;
			}
		}
		assert(j > -1);
		return static_cast<int>(j);
	}

	// Spherical Fibonacci Mapping by Benjamin Keinert et al. section 5.2.2
	static inline void SphFiboIndicesFromWorld(std::array<FiboInterpolateInfo, 9> * result, const Vec3 & p, const int nSamples)
	{
		assert(nSamples > 0);

		const auto madfrac = [](const double a, const double b) { return a * b - std::floor(a * b); };

		const double invNSamples = 1.0 / static_cast<double>(nSamples);
		const double GoldenPHI = (std::sqrt(5.0) + 1.0) / 2.0;
		const double theta = std::min(std::atan2(p[2], p[0]), Math::Pi);
		const double cosPhi = p[1];
		const double k = std::max(2.0, std::floor(std::log(nSamples * Math::Pi * std::sqrt(5.0) * (1.0 - cosPhi * cosPhi)) / std::log(GoldenPHI * GoldenPHI)));
		const double Fk = std::pow(GoldenPHI, k) / std::sqrt(5.0);
		const double F0 = std::round(Fk);
		const double F1 = std::round(Fk * GoldenPHI);
		const Vec2 BRow0(2.0 * Math::Pi * madfrac(F0 + 1.0, GoldenPHI - 1.0) - 2.0 * Math::Pi * (GoldenPHI - 1.0),
					  2.0 * Math::Pi * madfrac(F1 + 1.0, GoldenPHI - 1.0) - 2.0 * Math::Pi * (GoldenPHI - 1.0));
		const Vec2 BRow1(-2.0 * F0 * invNSamples,
						 -2.0 * F1 * invNSamples);
		const double invDetB = 1.0 / (BRow0[0] * BRow1[1] - BRow0[1] * BRow1[0]);
		const Vec2 invBRow0(BRow1[1] * invDetB, -BRow0[1] * invDetB);
		const Vec2 invBRow1(-BRow1[0] * invDetB, BRow0[0] * invDetB);
		const Vec2 cc(theta, cosPhi - (1.0 - 1.0 * invNSamples));
		const Vec2 cp = Math::Round(Vec2(Math::Dot(cc, invBRow0), Math::Dot(cc, invBRow1)));

		const int nInterp = 9;
		const int sqrtNInterp = 3;

		std::array<double, 9> distances;
		for (int s = 0; s < nInterp; s++)
		{
			const Vec2 uv(static_cast<double>(s % sqrtNInterp) - 1.0, static_cast<double>(s / sqrtNInterp) - 1.0);
			double qCosPhi = Math::Dot(BRow1, uv + cp) + (1.0 - invNSamples);
			qCosPhi = Math::Clamp(qCosPhi, -1.0, +1.0) * 2.0 - qCosPhi;

			const double i = std::floor(nSamples * 0.5 - qCosPhi * nSamples * 0.5);
			const double qTheta = 2.0 * Math::Pi * madfrac(i, GoldenPHI - 1.0);
			qCosPhi = 1.0 - (2.0 * i + 1) * invNSamples;
			const double qSinPhi = std::sqrt(1 - qCosPhi * qCosPhi);

			const Vec3 q(std::cos(qTheta) * qSinPhi, qCosPhi, std::sin(qTheta) * qSinPhi);
			const double squaredDistance = Math::Distance2(q, p);
			const double distance = std::sqrt(squaredDistance);
			distances[s] = distance;
			(*result)[s].index = static_cast<int>(i);
			assert((*result)[s].index >= 0);
		}

		// compute weight for each point
		const double h = std::sqrt(4.0 * Math::Pi * invNSamples);
		double sumWeight = 0.0;
		for (int s = 0; s < nInterp; s++)
		{
			const double ds = distances[s];
			const double t = std::max(0.0, 1.0 - ds / h);
			const double t2 = t * t;
			const double ws = 3.0*t2 - 2.0*t2*t;
			(*result)[s].weight = ws;
			sumWeight += ws;
		}

		// normalize all weight
		for (int s = 0; s < nInterp; s++) { (*result)[s].weight /= sumWeight; }
	}

	static inline Vec3 WorldFromPanorama(const Vec2 & uv)
	{
		double theta = uv[0] * Math::Pi * 2.0;
		double phi = uv[1] * Math::Pi;
		return WorldFromSpherical(Vec2(theta, phi));
	}

	static inline Vec2 World2FromAngle(const double theta)
	{
		return Vec2(std::sin(theta), std::cos(theta));
	}

	static inline Vec3 WorldFromSpherical(const Vec2 & thetaPhi)
	{
		const double & theta = thetaPhi[0];
		const double & phi = thetaPhi[1];
		const double sinphi = std::sin(phi);
		const double cosphi = std::cos(phi);
		const double sintheta = std::sin(theta);
		const double costheta = std::cos(theta);
		return Vec3(costheta * sinphi, cosphi, sintheta * sinphi);
	}

	static inline Vec3 WorldFromSphFiboIndex(const int index, const int nSamples)
	{
		assert(index >= 0);
		assert(nSamples >= 0);
		const double GoldenPHI = (std::sqrt(5.0) + 1.0) / 2.0;
		const double p = double(index) / GoldenPHI;
		const double q = p - std::floor(p);
		const double theta = 2.0 * Math::Pi * q;
		const double z = 1.0 - (2.0 * double(index) + 1.0) / double(nSamples);
		const double phi = std::acos(z);
		return WorldFromSpherical(Vec2(theta, phi));
	}
};