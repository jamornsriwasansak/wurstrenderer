#pragma once

#include "common/wurst.h"

#include "common/cdftable.h"
#include "common/geometry.h"
#include "common/texture.h"
#include "geometry/trianglemesh.h"

struct PhongMeshGeometry : public Geometry
{
	static void TwoPlanesFromRay(double * o1, Vec3 * d1, double * o2, Vec3 * d2, const Ray3 & ray)
	{
		CoordFrame3 cf(ray.direction, false);
		
		// first plane
		*d1 = Math::Cross(cf.x, cf.y);
		*o1 = -Math::Dot(ray.origin, *d1);

		// second plane
		*d2 = Math::Cross(cf.z, cf.y);
		*o2 = -Math::Dot(ray.origin, *d2);
	}

	static Vec3 EvalPhongTessellationPosition(const Vec3 & p0, const Vec3 & p1, const Vec3 & p2, const Vec3 & n0, const Vec3 & n1, const Vec3 & n2, const double u, const double v, const double alpha)
	{
		assert(u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && u + v <= 1.0);
		const Vec3 puv = u * p0 + v * p1 + (1.0 - u - v) * p2;
		const Vec3 pi0 = puv - Math::Dot(puv - p0, n0) * n0;
		const Vec3 pi1 = puv - Math::Dot(puv - p1, n1) * n1;
		const Vec3 pi2 = puv - Math::Dot(puv - p2, n2) * n2;
		const Vec3 pStar = u * pi0 + v * pi1 + (1.0 - u - v) * pi2;
		return (1.0 - alpha) * puv + alpha * pStar;
	}

	static Vec3 EvalPhongTessellationNormal(const Vec3 & n0, const Vec3 & n1, const Vec3 & n2, const double u, const double v, const Vec3 & c31, const Vec3 & c12, const Vec3 & c23, const Vec3 & e31, const Vec3 & e23)
	{
		// eq 29.
		const double w = 1.0 - u - v;
		const Vec3 dpdu = (w - u) * c31 + v * (c12 - c23) + e31;
		const Vec3 dpdv = (w - v) * c23 + u * (c12 - c31) - e23;
		const Vec3 n = Math::Cross(dpdu, dpdv);
		return Math::Normalize(n);
	}

	PhongMeshGeometry(const double alpha = 3.0 / 4.0):
		mPhongTessAlpha(alpha)
	{
	}

	// TODO:: implement
	double area() const override
	{
		return 0.0;
	}

	Bound3 bound() const override
	{
		Bound3 result;
		for (int i = 0; i < mPositions.size(); i++)
		{
			result = Bound3::Union(result, mPositions[i]);
		}
		return result;
	}

	Bound3 bound(const int primitiveIndex) const override
	{
		// fetch values and compute geometry normal
		Vec3 positions[3];
		Vec3 normals[3];
		for (int i = 0; i < 3; i++)
		{
			positions[i] = mPositions[mTriangles[primitiveIndex][i]];
			normals[i] = mNormals[mTriangles[primitiveIndex][i]];
		}
		Vec3 ng = Math::Normalize(Math::Cross(positions[1] - positions[0], positions[2] - positions[0]));

		// compute e and c values
		auto p_ = [&](const int i) -> Vec3& { return positions[i - 1]; };
		auto n_ = [&](const int i) -> Vec3& { return normals[i - 1]; };
		auto e_ = [&](int i, int j) -> Vec3 { return p_(j) - p_(i); };
		auto c_ = [&](int i, int j, Vec3 eij) -> Vec3 { return mPhongTessAlpha * (Math::Dot(n_(j), eij) * n_(j) - Math::Dot(n_(i), eij) * n_(i)); };
		const Vec3 e12 = e_(1, 2);
		const Vec3 e23 = e_(2, 3);
		const Vec3 e31 = e_(3, 1);
		const Vec3 c12 = c_(1, 2, e12);
		const Vec3 c23 = c_(2, 3, e23);
		const Vec3 c31 = c_(3, 1, e31);

		// compute u^ and v^ from eq 21, 22
		const Vec3 c31e31 = c31 + e31;
		const Vec3 c23e23 = c23 - e23;
		const Vec3 c122331 = c12 - c23 - c31;
		const double k = 1.0 / Math::Square(4.0 * Math::Dot(ng, c23) * Math::Dot(ng, c31) - Math::Dot(ng, c122331));
		const double uhat = k * (2.0 * Math::Dot(ng, c23) * Math::Dot(ng, c31e31) + Math::Dot(ng, c23e23) * Math::Dot(ng, c122331));
		const double vhat = k * (2.0 * Math::Dot(ng, c31) * Math::Dot(ng, c23e23) + Math::Dot(ng, c31e31) * Math::Dot(ng, c122331));

		// for computing distance from a point on implicit surface to the flat triangle
		auto d_ = [&](const double u, const double v) -> double
		{
			const Vec3 palpha = EvalPhongTessellationPosition(positions[0], positions[1], positions[2], normals[0], normals[1], normals[2], u, v, mPhongTessAlpha);
			return Math::Dot(ng, palpha - p_(1));
		};

		// depends on the case
		double thickness = 0.0;
		bool isUhatOkay = uhat >= 0.0 && uhat <= 1.0;
		bool isVhatOkay = vhat >= 0.0 && vhat <= 1.0;
		if (isUhatOkay && isVhatOkay && uhat + vhat <= 1.0)
		{
			thickness = d_(uhat, vhat);
		}
		else
		{
			double a = isUhatOkay ? d_(uhat, 0.0) : 0.0;
			double b = isVhatOkay ? d_(0.0, vhat) : 0.0;
			thickness = std::max(a, b);
		}
		assert(std::isfinite(thickness));

		// check whether the thickness is zero or not
		// if it is, simply return the triangle bound
		if (thickness == 0.0)
		{
			Bound3 triangleBound;
			for (int i = 0; i < 3; i++)
			{
				triangleBound = Bound3::Union(triangleBound, positions[i]);
			}
			return triangleBound;
		}

		// eq. 23
		auto sij_ = [&](const int i, const int j, const Vec3 & palphauv) -> double
		{
			return Math::Length(Math::Cross(p_(j) - palphauv, p_(i) - palphauv)) / Math::Length(e_(i, j));
		};

		// for all edges
		double sidedrops[3];
		sidedrops[0] = 0.0; sidedrops[1] = 0.0; sidedrops[2] = 0.0;
		const std::array<Vec3, 3> uvws = { Vec3(0.0, 0.5, 0.5), Vec3(0.5, 0.0, 0.5), Vec3(0.5, 0.5, 0.0) };
		for (int iMidpoint = 0; iMidpoint < 3; iMidpoint++)
		{
			const Vec3 & uvw = uvws[iMidpoint];
			const Vec3 palphauv = EvalPhongTessellationPosition(positions[0], positions[1], positions[2], normals[0], normals[1], normals[2], uvw[0], uvw[1], mPhongTessAlpha);
			// compute sidedrop
			for (int iEdge = 0; iEdge < 3; iEdge++)
			{
				const double slength = sij_(iEdge + 1, (iEdge + 1) % 3 + 1, palphauv);
				sidedrops[iEdge] = std::max(slength, sidedrops[iEdge]);
				assert(std::isfinite(sidedrops[iEdge]));
			}
		}

		Bound3 result;
		for (int iEdge = 0; iEdge < 3; iEdge++)
		{
			const int i = iEdge;
			const int j = (iEdge + 1) % 3;

			// move the edge by the distance of sidedrop and merge into result
			const Vec3 eij = positions[j] - positions[i];

			// build orthonormal basis 
			const Vec3 & y = ng;
			const Vec3 x = Math::Normalize(eij);
			const Vec3 z = Math::Cross(x, y);

			// move both vertices on end of this edge toward -z
			const Vec3 vi = z * sidedrops[i] + positions[i];
			const Vec3 vj = z * sidedrops[i] + positions[j];

			result = Bound3::Union(result, vi);
			result = Bound3::Union(result, vj);
			result = Bound3::Union(result, vi + y * thickness);
			result = Bound3::Union(result, vj + y * thickness);
		}
		return result;
	}

	double evalPdfA() const override
	{
		throw std::runtime_error("unimpl");
		return 0.0;
	}

	int getNumPrimitives() const override
	{
		return int(mTriangles.size());
	}

	bool intersect(SimpleVertex * isect, Ray3 * ray, const int primitiveIndex) const override
	{
		int axis = int(Math::MaxDimension(ray->direction));

		// fetch values
		Vec3 positions[3];
		Vec3 normals[3];
		for (int i = 0; i < 3; i++)
		{
			positions[i] = mPositions[mTriangles[primitiveIndex][i]];
			normals[i] = mNormals[mTriangles[primitiveIndex][i]];
		}
		auto p_ = [&](const int i) -> Vec3& { return positions[i - 1]; };
		auto n_ = [&](const int i) -> Vec3& { return normals[i - 1]; };

		// TODO: choose a plane better so that a = 0 and m = 0
		// decompose the ray into the intersection of two planes
		double o1, o2;
		Vec3 d1, d2;
		TwoPlanesFromRay(&o1, &d1, &o2, &d2, *ray);

		// solve for a b c d e f l m n o p q - equation 5, 6
		auto e_ = [&](int i, int j) -> Vec3 { return p_(j) - p_(i); };
		auto c_ = [&](int i, int j, Vec3 eij) -> Vec3 { return mPhongTessAlpha * (Math::Dot(n_(j), eij) * n_(j) - Math::Dot(n_(i), eij) * n_(i)); };

		const Vec3 e12 = e_(1, 2);
		const Vec3 e23 = e_(2, 3);
		const Vec3 e31 = e_(3, 1);

		const Vec3 c12 = c_(1, 2, e12);
		const Vec3 c23 = c_(2, 3, e23);
		const Vec3 c31 = c_(3, 1, e31);

		const Vec3 c123 = c12 - c23 - c31;
		const Vec3 ce31 = c31 + e31;
		const Vec3 ce23 = c23 - e23;

		const double a = -Math::Dot(d1, c31);
		const double b = -Math::Dot(d1, c23);
		const double c = Math::Dot(d1, p_(3)) + o1;
		const double d = Math::Dot(d1, c123);
		const double e = Math::Dot(d1, ce31);
		const double f = Math::Dot(d1, ce23);
		const double l = -Math::Dot(d2, c31);
		const double m = -Math::Dot(d2, c23);
		const double n = Math::Dot(d2, p_(3)) + o2;
		const double o = Math::Dot(d2, c123);
		const double p = Math::Dot(d2, ce31);
		const double q = Math::Dot(d2, ce23);

		// solve for a3 a2 a1 a0 - equation 16
		const double a3 = a * b*c + (d*e*f - a * f*f - b * e*e - c * d*d) * 0.25;
		const double a2 = a * b*n + a * m*c + l * b*c - (a*f*q + b * e*p + c * d*o) * 0.5 + (o*e*f + d * e*q + d * p*f - l * f*f - m * e*e - n * d*d) * 0.25;
		const double a1 = a * m*n + l * b*n + l * m*c - (l*f*q + m * e*p + n * d*o) * 0.5 + (d*p*q + o * e*q + o * p*f - a * q*q - b * p*p - c * o*o) * 0.25;
		const double a0 = l * m*n + (o*p*q - l * q*q - m * p*p - n * o*o) * 0.25;

		// solve for x
		double xs[3];
		int numXs = Math::SolveRealCubic(&xs[0], &xs[1], &xs[2], a3, a2, a1, a0);
		if (numXs == 0) return false;

		// (based on intersection test subsection) choose best x that minimize minXCost = M12*M21 - M11*M22 = M12*M12 - M11*M22
		double x;
		double minXCost = std::numeric_limits<double>::infinity();
		for (int i = 0; i < numXs; i++)
		{
			const double m11 = xs[i] * a + l;
			const double m22 = xs[i] * b + m;
			const double m12 = (xs[i] * d + o) * 0.5;
			const double xCost = m12 * m12 - m11 * m22;
			if (xCost < minXCost)
			{
				x = xs[i];
				minXCost = xCost;
			}
		}

		// (at the end of section 3.2) terminate if minXCost = M12&M21 - M11*M22 < 0
		if (minXCost < 0) return false;

		// determine the case for more robustness
		const double m11 = x * a + l;
		const double m22 = x * b + m;
		const bool isGreaterCase = std::abs(m11) > std::abs(m22);

		// l1 = alpha_1*u  + beta_1*v + gamma_1 = 0
		// l2 = alpha_2*u  + beta_2*v + gamma_2 = 0
		// since either alpha or beta is 1.0 depends on the case, we can save some space by sharing memory between alpha and beta
		// we call it albeta
		double albeta[2]; 
		double gamma[2];
		if (isGreaterCase) // abs(m11) > abs(m22), in the paper they forgot abs
		{
			// alpha is 1.0
			// this albeta is beta
			const double invM11 = 1.0 / m11;
			const double mp22 = m22 * invM11;
			const double mp33 = (x * c + n) * invM11;
			const double mp12 = (x * d + o) * 0.5 * invM11;
			const double mp13 = (x * e + p) * 0.5 * invM11;
			const double mp23 = (x * f + q) * 0.5 * invM11;
			const double betaSqrtTerm = std::sqrt(mp12 * mp12 - mp22);
			const double gammaSqrtTerm = std::sqrt(mp13 * mp13 - mp33);
			albeta[0] = mp12 + betaSqrtTerm;
			albeta[1] = mp12 - betaSqrtTerm;
			gamma[0] = mp13 + gammaSqrtTerm;
			gamma[1] = mp13 - gammaSqrtTerm;

			const double diffA = std::abs((2.0 * mp23) - (albeta[0] * gamma[0] + albeta[1] * gamma[1]));
			const double diffB = std::abs((2.0 * mp23) - (albeta[0] * gamma[1] + albeta[1] * gamma[0]));
			if (diffA < diffB) { std::swap(gamma[0], gamma[1]); }
		}
		else // abs(m22) >= abs(m11)
		{
			// beta is 1.0
			// this albeta is alpha
			const double invM22 = 1.0 / m22;
			const double mp11 = m11 * invM22;
			const double mp12 = (x * d + o) * 0.5 * invM22;
			const double mp13 = (x * e + p) * 0.5 * invM22;
			const double mp23 = (x * f + q) * 0.5 * invM22;
			const double mp33 = (x * c + n) * invM22;
			const double alphaSqrtTerm = std::sqrt(mp12 * mp12 - mp11);
			const double gammaSqrtTerm = std::sqrt(mp23 * mp23 - mp33);
			albeta[0] = mp12 + alphaSqrtTerm;
			albeta[1] = mp12 - alphaSqrtTerm;
			gamma[0] = mp23 + gammaSqrtTerm;
			gamma[1] = mp23 - gammaSqrtTerm;

			const double diffA = std::abs((2.0 * mp13) - (albeta[0] * gamma[0] + albeta[1] * gamma[1]));
			const double diffB = std::abs((2.0 * mp13) - (albeta[0] * gamma[1] + albeta[1] * gamma[0]));
			if (diffA < diffB) { std::swap(gamma[0], gamma[1]); }
		}

		// since l1 * l2 = 0 from eq. 14 = eq. 10 there are two possible cases. l1 = 0 or l2 = 0 we thus have to solve twice
		bool didIntersect = false;
		for (int i = 0; i < 2; i++)
		{
			const double g = -albeta[i];
			const double h = -gamma[i];

			double k2, k1, k0;
			if (isGreaterCase)
			{
				// u = gv + h = -beta*v - gamma
				// substitute u into eq 5., we then have a quadratic equation
				k2 = a*g*g + b + d*g;
				k1 = 2.0*a*g*h + d*h + e*g + f;
				k0 = a*h*h + c + e*h;
			}
			else
			{
				// v = gu + h = -alpha*v - gamma
				// substitute v into eq 5., we then have a quadratic equation
				k2 = a + b*g*g + d*g;
				k1 = 2.0*b*g*h + d*h + e + f*g;
				k0 = b*h*h + c + f*h;
			}

			double ans[2];
			if (Math::SolveQuadratic(&ans[0], &ans[1], k2, k1, k0))
			{
				for (int iAns = 0; iAns < 2; iAns++)
				{
					double u, v;
					if (isGreaterCase)
					{
						v = ans[iAns];
						u = g * v + h;
					}
					else
					{
						u = ans[iAns];
						v = g * u + h;
					}
					const double w = 1.0 - u - v;
					if (v >= 0.0 && u >= 0.0 && w >= 0.0)
					{
						if (isect)
						{
							const Vec3 position = EvalPhongTessellationPosition(positions[0], positions[1], positions[2], normals[0], normals[1], normals[2], u, v, mPhongTessAlpha);
							const double t = (position[axis] - ray->origin[axis]) / ray->direction[axis];
							if (ray->tMax >= t && ray->tMin <= t)
							{
								didIntersect = true;
								ray->tMax = t;
								isect->mGeometryNormal = EvalPhongTessellationNormal(normals[0], normals[1], normals[2], u, v, c31, c12, c23, e31, e23);
							}
						}
					}
				}
			}
		}

		return didIntersect;
	}

	void precomputeForSampling() override
	{
		throw std::runtime_error("unimpl");
	}

	SimpleVertex sample(double * pdfA, const Vec2 & sample) const override
	{
		throw std::runtime_error("unimpl");
		return SimpleVertex();
	}

	void transform(const Matrix4 & matrix) override
	{
		throw std::runtime_error("unimpl");
	}

	double mPhongTessAlpha;

	// arguments for 3d

	std::vector<Vec3> mPositions;
	std::vector<Vec3> mNormals;
	std::vector<Vec2> mTextureCoords;
	std::vector<Ivec3> mTriangles;
};
