#pragma once

#include "scalarmath.h"

#include <functional>
#include <cassert>
#include <sstream>

/*

// else always true
inline bool isfinite2(int x) { return true; }
inline bool isfinite2(unsigned int) { return true; }
inline bool isfinite2(size_t) { return true; }
*/

namespace Math
{
	template <const size_t Size, const size_t Id>
	struct alignas(32) VecSimd4;

	template <typename Scalar, const size_t Size, const size_t Id>
	struct VecPrim
	{
		static const size_t NumElements = Size;

		// printing stuffs

		friend std::ostream & operator<<(std::ostream & out, const VecPrim & v)
		{
			out << "(";
			for (size_t i = 0; i < Size - 1; i++) { out << v[i] << ", "; }
			out << v[Size - 1] << ")";
			return out;
		}

		/*
		template <typename Scalar, const size_t Size, const size_t Id>
		friend std::ostream & operator<<(std::ostream & out, const std::vector<VecPrim<Scalar, Size, Id>> & vs)
		{
			const Uint numDisplay = 6;
			out << "[";
			for (Uint i = 0; i < vs.size(); i++)
			{
				if (i < Math::DivAndCeil(numDisplay, 2) || i + numDisplay / 2 >= vs.size())
				{
					out << vs[i];
					if (i != vs.size() - 1) out << " ";
				}
				else if (vs.size() > numDisplay && i == Math::DivAndCeil(numDisplay, 2))
				{
					out << "...";
					if (i != vs.size() - 1) out << " ";
				}
			}
			out << "]";
			return out;
		}
		*/

		// conversion VecSimd4 -> VecPrim

		template <size_t OtherSize, size_t OtherId>
		__forceinline explicit VecPrim(const VecSimd4<OtherSize, OtherId> & u);

		// VecPrim implementation

		__forceinline VecPrim()
		{
			for (size_t i = 0; i < Size; i++) v[i] = 0;
			assert(IsFinite(*this));
		}

		__forceinline explicit VecPrim(const Scalar s)
		{
			for (size_t i = 0; i < Size; i++) v[i] = s;
			assert(IsFinite(*this));
		}

		__forceinline VecPrim(const VecPrim & u)
		{
			for (size_t i = 0; i < Size; i++) v[i] = u[i];
			assert(IsFinite(*this));
		}

		template <typename OtherScalar, size_t OtherSize, size_t OtherId>
		__forceinline explicit VecPrim(const VecPrim<OtherScalar, OtherSize, OtherId> & u)
		{
			static_assert(OtherSize == Size, "VecPrim: Must initialize with the vector of equal Size");
			for (size_t i = 0; i < OtherSize; i++) v[i] = static_cast<Scalar>(u[i]);
			assert(IsFinite(*this));
		}

		template <typename ... MoreScalars>
		__forceinline explicit VecPrim(const Scalar s, MoreScalars ... t)
		{
			static_assert(sizeof...(MoreScalars) == Size - 1, "VecPrim: number of arguments mismatch the Size");
			variadicInitialize(0, s, t...);
			assert(IsFinite(*this));
		}

		template <typename ... MoreScalars>
		__forceinline void variadicInitialize(const size_t index, const Scalar s, MoreScalars ... t)
		{
			//static_assert(std::is_floating_point<Scalar>::value || std::is_arithmetic<Scalar>::value);
			v[index] = s;
			variadicInitialize(index + 1, t...);
		}

		__forceinline void variadicInitialize(const size_t index, const Scalar s)
		{
			v[index] = s;
		}

		__forceinline const Scalar operator[](const size_t index) const
		{
			assert(index < Size);
			return v[index];
		}

		__forceinline Scalar & operator[](const size_t index)
		{
			assert(index < Size);
			return v[index];
		}

		__forceinline VecPrim operate(const VecPrim & other, const std::function<Scalar(const Scalar a, const Scalar b)>& func) const
		{
			VecPrim result;
			for (size_t i = 0; i < Size; i++) { result[i] = func(v[i], other[i]); }
			assert(IsFinite(*this));
			return result;
		}

		__forceinline VecPrim operator+(const VecPrim & other) const { return operate(other, [](const Scalar a, const Scalar b) { return a + b; }); }
		__forceinline VecPrim operator-(const VecPrim & other) const { return operate(other, [](const Scalar a, const Scalar b) { return a - b; }); }
		__forceinline VecPrim operator*(const VecPrim & other) const { return operate(other, [](const Scalar a, const Scalar b) { return a * b; }); }
		__forceinline VecPrim operator/(const VecPrim & other) const { return operate(other, [](const Scalar a, const Scalar b) { return a / b; }); }

		__forceinline VecPrim operate(const std::function<Scalar(const Scalar a)>& func) const
		{
			VecPrim result;
			for (size_t i = 0; i < Size; i++) { result[i] = func(v[i]); }
			assert(IsFinite(*this));
			return result;
		}

		__forceinline VecPrim operator-() const { return operate([](const Scalar a) { return -a; }); }

		__forceinline VecPrim operator+(const Scalar other) const { return operate([&](const Scalar a) { return a + other; }); }
		__forceinline VecPrim operator-(const Scalar other) const { return operate([&](const Scalar a) { return a - other; }); }
		__forceinline VecPrim operator*(const Scalar other) const { return operate([&](const Scalar a) { return a * other; }); }
		__forceinline VecPrim operator/(const Scalar other) const { return operate([&](const Scalar a) { return a / other; }); }

		__forceinline friend VecPrim operator+(const Scalar other, const VecPrim & v) { return v + other; }
		__forceinline friend VecPrim operator-(const Scalar other, const VecPrim & v) { return v - other; }
		__forceinline friend VecPrim operator*(const Scalar other, const VecPrim & v) { return v * other; }
		__forceinline friend VecPrim operator/(const Scalar other, const VecPrim & v) { return v.operate([&](const Scalar a) { return other / a; }); }

		__forceinline VecPrim & selfOperate(const VecPrim & other, const std::function<void(Scalar & a, const Scalar & b)>& func)
		{
			for (size_t i = 0; i < Size; i++) { func(v[i], other[i]); }
			assert(IsFinite(*this));
			return *this;
		}

		__forceinline VecPrim & operator+=(const VecPrim & other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a += b; }); }
		__forceinline VecPrim & operator-=(const VecPrim & other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a -= b; }); }
		__forceinline VecPrim & operator*=(const VecPrim & other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a *= b; }); }
		__forceinline VecPrim & operator/=(const VecPrim & other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a /= b; }); }

		__forceinline VecPrim & selfOperate(const Scalar & other, const std::function<void(Scalar & a, const Scalar & b)>& func)
		{
			for (size_t i = 0; i < Size; i++) { func(v[i], other); }
			assert(IsFinite(*this));
			return *this;
		}

		__forceinline VecPrim & operator+=(const Scalar other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a += b; }); }
		__forceinline VecPrim & operator-=(const Scalar other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a -= b; }); }
		__forceinline VecPrim & operator*=(const Scalar other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a *= b; }); }
		__forceinline VecPrim & operator/=(const Scalar other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a /= b; }); }

		__forceinline bool operator==(const VecPrim & other) const { for (size_t i = 0; i < Size; i++) { if (v[i] != other[i]) return false; } return true; }
		__forceinline bool operator!=(const VecPrim & other) const { return !(*this == other); }

		Scalar v[Size];
	};

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Abs(const VecPrim<Scalar, Size, Id> & v)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::abs(v[i]);
		assert(IsFinite(result));
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline Scalar AbsDot(const VecPrim<Scalar, Size, Id> & u, const VecPrim<Scalar, Size, Id> & v)
	{
		Scalar result = static_cast<Scalar>(0);
		for (size_t i = 0; i < Size; i++) result += u[i] * v[i];
		return std::abs(result);
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline bool AnyGreater(const VecPrim<Scalar, Size, Id> & a, const VecPrim<Scalar, Size, Id> & b)
	{
		for (size_t i = 0; i < Size; i++) if (a[i] > b[i]) return true;
		return false;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline Scalar ClampDot(const VecPrim<Scalar, Size, Id> & u, const VecPrim<Scalar, Size, Id> & v)
	{
		Scalar result = static_cast<Scalar>(0);
		for (size_t i = 0; i < Size; i++) result += u[i] * v[i];
		return std::max(result, static_cast<Scalar>(0));
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Ceil(const VecPrim<Scalar, Size, Id> & v)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::ceil(v[i]);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Clamp(const VecPrim<Scalar, Size, Id> & v, const VecPrim<Scalar, Size, Id> & a, const VecPrim<Scalar, Size, Id> & b)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::min(std::max(v[i], a[i]), b[i]);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Cross(const VecPrim<Scalar, Size, Id> & u, const VecPrim<Scalar, Size, Id> & v)
	{
		return VecPrim<Scalar, Size, Id>(u[1] * v[2] - v[1] * u[2], u[2] * v[0] - v[2] * u[0], u[0] * v[1] - v[0] * u[1]);
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline Scalar Distance(const VecPrim<Scalar, Size, Id> & u, const VecPrim<Scalar, Size, Id> & v)
	{
		return std::sqrt(Distance2(u, v));
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline Scalar Distance2(const VecPrim<Scalar, Size, Id> & u, const VecPrim<Scalar, Size, Id> & v)
	{
		return Length2(u - v);
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline Scalar Dot(const VecPrim<Scalar, Size, Id> & u, const VecPrim<Scalar, Size, Id> & v)
	{
		Scalar result = 0;
		for (size_t i = 0; i < Size; i++) result += u[i] * v[i];
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Exp(const VecPrim<Scalar, Size, Id> & v)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::exp(v[i]);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Floor(const VecPrim<Scalar, Size, Id> & v)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::floor(v[i]);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline bool IsNormalized(const VecPrim<Scalar, Size, Id> & v)
	{
		return Math::IsApprox(Length2(v), 1.0);
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline bool IsZero(const VecPrim<Scalar, Size, Id> & v)
	{
		for (size_t i = 0; i < Size; i++) if (v[i] != 0.0) return false;
		return true;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline bool IsFinite(const VecPrim<Scalar, Size, Id> & v)
	{
		for (size_t i = 0; i < Size; i++) if (!isfinite2(v[i])) return false;
		return true;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline Scalar Length(const VecPrim<Scalar, Size, Id> & v)
	{
		return std::sqrt(Length2(v));
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline Scalar Length2(const VecPrim<Scalar, Size, Id> & v)
	{
		return Dot(v, v);
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Pow(const VecPrim<Scalar, Size, Id> & v, const Scalar & s)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::pow(v[i], s);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Pow(const VecPrim<Scalar, Size, Id> & u, const VecPrim<Scalar, Size, Id> & v)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::pow(u[i], v[i]);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline size_t MaxDimension(const VecPrim<Scalar, Size, Id> & v)
	{
		size_t index = 0;
		Scalar max = v[0];
		for (size_t i = 1; i < Size; i++) if (v[i] > max) { index = i; max = v[i]; }
		return index;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Max(const VecPrim<Scalar, Size, Id> & u, const VecPrim<Scalar, Size, Id> & v)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::max(u[i], v[i]);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Min(const VecPrim<Scalar, Size, Id> & u, const VecPrim<Scalar, Size, Id> & v)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::min(u[i], v[i]);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Normalize(const VecPrim<Scalar, Size, Id> & v)
	{
		return v / Length(v);
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Round(const VecPrim<Scalar, Size, Id> & v)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::round(v[i]);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline VecPrim<Scalar, Size, Id> Sqrt(const VecPrim<Scalar, Size, Id> & v)
	{
		VecPrim<Scalar, Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::sqrt(v[i]);
		return result;
	}

	template <typename Scalar, const size_t Size, const size_t Id>
	__forceinline Scalar Volume(const VecPrim<Scalar, Size, Id> & v)
	{
		Scalar result = static_cast<Scalar>(1);
		for (size_t i = 0; i < Size; i++) result *= v[i];
		return result;
	}
};