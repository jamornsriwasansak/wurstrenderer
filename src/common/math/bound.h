#pragma once

#include <numeric>
#include "vecmath.h"

template <typename Vec, typename Scalar>
struct Bound
{
	friend std::ostream & operator<<(std::ostream & out, const Bound & v)
	{
		out << "[" << v.pMin << ", " << v.pMax << "]";
		return out;
	}

	static Scalar Area(const Bound & a)
	{
		Scalar result(0);
		Vec size = Math::Max(a.pMax - a.pMin, Vec3(0.0));
		for (Uint i = 0; i < Vec::NumElements; i++)
			for (Uint j = i + 1; j < Vec::NumElements; j++)
			{
				result += size[i] * size[j] * Scalar(Vec::NumElements - 1);
			}
		return result;
	}

	static Bound Intersect(const Bound & a, const Bound & b)
	{
		return Bound(Math::Max(a.pMin, b.pMin), Math::Min(a.pMax, b.pMax));
	}

	static bool IsInside(const Bound & a, const Vec & b)
	{
		for (Uint i = 0; i < Vec::NumElements; i++)
		{
			if (b[i] < a.pMin[i] || b[i] >= a.pMax[i]) { return false; }
		}
		return true;
	}

	static Bound Padding(const Bound & a, const double padSize)
	{
		return Bound(a.pMin - Vec(padSize), a.pMax + Vec(padSize));
	}

	static Bound Resize(const Bound & a, const double ratio)
	{
		assert(!(a.pMin == Vec(std::numeric_limits<Scalar>::max())));
		assert(!(a.pMax == Vec(std::numeric_limits<Scalar>::lowest())));
		Vec center = (a.pMin + a.pMax) * 0.5;
		Vec halfSize = (a.pMax - a.pMin) * ratio * 0.5;
		return Bound(center - halfSize, center + halfSize);
	}

	static Bound Union(const Bound & a, const Vec & b)
	{
		return Bound(Math::Min(a.pMin, b), Math::Max(a.pMax, b));
	}

	static Bound Union(const Bound & a, const Bound & b)
	{
		return Bound(Math::Min(a.pMin, b.pMin), Math::Max(a.pMax, b.pMax));
	}

	static int Volume(const Bound & a)
	{
		return Math::Volume(Math::Max(a.pMax - a.pMin, Ivec2(0)));
	}

	Bound() : pMin(std::numeric_limits<double>::max()), pMax(std::numeric_limits<double>::lowest()) {}

	Bound(const Vec & p) : pMin(p), pMax(p)
	{
	}

	Bound(const Vec & min, const Vec & max): pMin(min), pMax(max)
	{
		//assert(!Math::AnyGreater(min, max));
	}

	Vec pMin = Vec(std::numeric_limits<Scalar>::max());
	Vec pMax = Vec(std::numeric_limits<Scalar>::lowest());
};

using Bound3 = Bound<Vec3, double>;
using Bound2 = Bound<Vec2, double>;
using Ibound3 = Bound<Ivec3, int>;
using Ibound2 = Bound<Ivec2, int>;
using Ubound3 = Bound<Uvec3, Uint>;
using Ubound2 = Bound<Uvec2, Uint>;

namespace Math
{
	static int Index(const Ibound2 & b, const Ivec2 & p)
	{
		assert(Ibound2::IsInside(b, p));
		const Ivec2 vol = b.pMax - b.pMin;
		const Ivec2 diff = p - b.pMin;
		return diff[1] * vol[0] + diff[0];
	}
};

template <typename Bound, typename Vec, typename Scalar>
class Bound2Iterator : public std::forward_iterator_tag {
public:
	Bound2Iterator(const Bound &b, const Vec &pt) : mCurrent(pt), mBound(&b)
	{
	}

	Bound2Iterator operator++()
	{
		advance();
		return *this;
	}

	Bound2Iterator operator++(int)
	{
		Bound2Iterator old = *this;
		advance();
		return old;
	}

	bool operator==(const Bound2Iterator &iter) const
	{
		return mCurrent == iter.mCurrent && mBound == iter.mBound;
	}

	bool operator!=(const Bound2Iterator &iter) const
	{
		return mCurrent != iter.mCurrent || mBound != iter.mBound;
	}

	Ivec2 operator*() const
	{
		return mCurrent;
	}

	void advance()
	{
		if (++mCurrent[0] == mBound->pMax[0])
		{
			mCurrent[0] = mBound->pMin[0];
			++mCurrent[1];
		}
	}

	Vec mCurrent;
	const Bound * mBound;
};

inline Bound2Iterator<Ibound2, Ivec2, int> begin(const Ibound2 & bound)
{
	return Bound2Iterator<Ibound2, Ivec2, int>(bound, bound.pMin);
}

inline Bound2Iterator<Ibound2, Ivec2, int> end(const Ibound2 & bound)
{
	Ivec2 pEnd(bound.pMin[0], bound.pMax[1]);
	if (bound.pMin[0] >= bound.pMax[0] || bound.pMin[1] >= bound.pMax[1]) { pEnd = bound.pMin; }
	return Bound2Iterator<Ibound2, Ivec2, int>(bound, pEnd);
}

inline Bound2Iterator<Ubound2, Uvec2, Uint> begin(const Ubound2 & bound)
{
	return Bound2Iterator<Ubound2, Uvec2, Uint>(bound, bound.pMin);
}

inline Bound2Iterator<Ubound2, Uvec2, Uint> end(const Ubound2 & bound)
{
	Uvec2 pEnd(bound.pMin[0], bound.pMax[1]);
	if (bound.pMin[0] >= bound.pMax[0] || bound.pMin[1] >= bound.pMax[1]) { pEnd = bound.pMin; }
	return Bound2Iterator<Ubound2, Uvec2, Uint>(bound, pEnd);
}
