#pragma once

#include "scalarmath.h"

#include <algorithm>
#include <cassert>
#include <immintrin.h>

namespace Math
{
	template <typename Scalar, const size_t Size, const size_t Id>
	struct VecPrim;

	template <const size_t Size, const size_t Id>
	struct alignas(32) VecSimd4
	{
		static const size_t NumElements = Size;

		// printing stuffs

		friend std::ostream & operator<<(std::ostream & out, const VecSimd4 & v)
		{
			out << "(";
			for (size_t i = 0; i < Size - 1; i++) { out << v[i] << ", "; }
			out << v[Size - 1] << ")";
			return out;
		}

		/*
		template <const size_t Size, const size_t Id>
		friend std::ostream & operator<<(std::ostream & out, const std::vector<VecSimd4<Size, Id>> & vs)
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

		// conversion VecPrim -> VecSimd4

		template <typename OtherScalar, const size_t OtherSize, const size_t OtherId>
		__forceinline explicit VecSimd4(const VecPrim<OtherScalar, OtherSize, OtherId> & ov);

		// VecSimd4 implementation

		__forceinline VecSimd4() : v(_mm256_setzero_pd()) {}

		__forceinline explicit VecSimd4(const double a) : v((Size == 3) ? _mm256_setr_pd(a, a, a, 0) : _mm256_set1_pd(a))
		{
			assert(IsFinite(*this));
		}

		__forceinline VecSimd4(const double a, const double b, const double c) : v(_mm256_setr_pd(a, b, c, 0))
		{
			static_assert(Size == 3, "VecSimd4: 3 arguments constructor is for Size of 3 only");
			assert(IsFinite(*this));
		}

		__forceinline VecSimd4(const double a, const double b, const double c, const double d) : v(_mm256_setr_pd(a, b, c, d))
		{
			static_assert(Size == 4, "VecSimd4: 4 arguments constructor is for Size of 4 only");
			assert(IsFinite(*this));
		}

		__forceinline VecSimd4(const VecSimd4 & ov) : v(ov.v)
		{
			assert(IsFinite(*this));
		}

		__forceinline VecSimd4(const __m256d & ov) : v(ov)
		{
			assert(IsFinite(*this));
		}

		template <const size_t OtherSize, const size_t OtherId>
		__forceinline explicit VecSimd4(const VecSimd4<OtherSize, OtherId> & ov) : v(ov.v)
		{
			static_assert(Size <= OtherSize, "VecSimd4: the information is not enough to create vec of this sizes");
			assert(IsFinite(*this));
		}

		const double operator[](const Uint index) const
		{
			assert(index < Size);
			return u[index];
		}

		double & operator[](const Uint index)
		{
			assert(index < Size);
			return u[index];
		}

		void * operator new(const size_t size)
		{
			return _mm_malloc(size, 32);
		}

		void operator delete(void * p)
		{
			return _mm_free(p);
		}

		friend __forceinline VecSimd4 operator+(const VecSimd4 & a, const VecSimd4 & b) { return _mm256_add_pd(a.v, b.v); }
		friend __forceinline VecSimd4 operator+(const VecSimd4 & a, const double b) { return a + VecSimd4(b); }
		friend __forceinline VecSimd4 operator+(const double a, const VecSimd4 & b) { return VecSimd4(a) + b; }

		friend __forceinline VecSimd4 operator-(const VecSimd4 & a, const VecSimd4 & b) { return _mm256_sub_pd(a.v, b.v); }
		friend __forceinline VecSimd4 operator-(const VecSimd4 & a, const double b) { return a - VecSimd4(b); }
		friend __forceinline VecSimd4 operator-(const double a, const VecSimd4 & b) { return VecSimd4(a) - b; }
		friend __forceinline VecSimd4 operator-(const VecSimd4 & a) { return _mm256_xor_pd(a.v, _mm256_set1_pd(-0.)); }

		friend __forceinline VecSimd4 operator*(const VecSimd4 & a, const VecSimd4 & b) { return _mm256_mul_pd(a.v, b.v); }
		friend __forceinline VecSimd4 operator*(const VecSimd4 & a, const double b) { return a * VecSimd4(b); }
		friend __forceinline VecSimd4 operator*(const double a, const VecSimd4 & b) { return VecSimd4(a) * b; }

		friend __forceinline VecSimd4 operator/(const VecSimd4 & a, const VecSimd4 & b) { return _mm256_div_pd(a.v, b.v); }
		friend __forceinline VecSimd4 operator/(const VecSimd4 & a, const double b) { return a / VecSimd4(b); }
		friend __forceinline VecSimd4 operator/(const double a, const VecSimd4 & b) { return VecSimd4(a) / b; }

		friend __forceinline VecSimd4 operator+=(VecSimd4 & a, const VecSimd4 & b) { return a = a + b; }
		friend __forceinline VecSimd4 operator+=(VecSimd4 & a, const double b) { return a = a + b; }

		friend __forceinline VecSimd4 operator-=(VecSimd4 & a, const VecSimd4 & b) { return a = a - b; }
		friend __forceinline VecSimd4 operator-=(VecSimd4 & a, const double b) { return a = a - b; }

		friend __forceinline VecSimd4 operator*=(VecSimd4 & a, const VecSimd4 & b) { return a = a * b; }
		friend __forceinline VecSimd4 operator*=(VecSimd4 & a, const double b) { return a = a * b; }

		friend __forceinline VecSimd4 operator/=(VecSimd4 & a, const VecSimd4 & b) { return a = a / b; }
		friend __forceinline VecSimd4 operator/=(VecSimd4 & a, const double b) { return a = a / b; }

		friend __forceinline bool operator==(const VecSimd4 & a, const VecSimd4 & b)
		{
			if (Size == 4)
			{
				__m256d cmp = _mm256_cmp_pd(a.v, b.v, _CMP_EQ_US);
				return _mm256_movemask_pd(cmp) == 0b1111;
			}
			else
			{
				__m256d cmp = _mm256_cmp_pd(a.v, b.v, _CMP_EQ_US);
				return (_mm256_movemask_pd(cmp) & 0b0111) == 0b0111;
			}
		}

		union
		{
			double u[4];
			__m256d v;
		};
	};

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Abs(const VecSimd4<Size, Id> & v)
	{
		return _mm256_andnot_pd(_mm256_set1_pd(-0.), v.v);
	}

	template <const size_t Size, const size_t Id>
	__forceinline double AbsDot(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		return std::abs(Dot(u, v));
	}

	template <const size_t Size, const size_t Id>
	__forceinline bool AnyGreater(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		for (size_t i = 0; i < Size; i++) if (u[i] > v[i]) return true;
		return false;
	}

	template <const size_t Size, const size_t Id>
	__forceinline double Dot(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		double sum = u[0] * v[0];
		for (Uint i = 1; i < Size; i++) { sum += u[i] * v[i]; }
		return sum;
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Ceil(const VecSimd4<Size, Id> & v)
	{
		return _mm256_ceil_pd(v.v);
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Clamp(const VecSimd4<Size, Id> & v, const VecSimd4<Size, Id> & a, const VecSimd4<Size, Id> & b)
	{
		return Max(Min(v, b), a);
	}

	template <const size_t Size, const size_t Id>
	__forceinline double ClampDot(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		return std::max(Dot(u, v), 0.0);
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Cross(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		static_assert(Size == 3, "cross is only allowed when size == 3 only");
		__m256d u120 = _mm256_permute4x64_pd(u.v, _MM_SHUFFLE(3, 0, 2, 1));
		__m256d v201 = _mm256_permute4x64_pd(v.v, _MM_SHUFFLE(3, 1, 0, 2));
		__m256d uv1 = _mm256_mul_pd(u120, v201);
		__m256d u201 = _mm256_permute4x64_pd(u.v, _MM_SHUFFLE(3, 1, 0, 2));
		__m256d v120 = _mm256_permute4x64_pd(v.v, _MM_SHUFFLE(3, 0, 2, 1));
		__m256d uv2 = _mm256_mul_pd(u201, v120);
		return _mm256_sub_pd(uv1, uv2);
	}

	template <const size_t Size, const size_t Id>
	__forceinline double Distance(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		return std::sqrt(Distance2(u, v));
	}

	template <const size_t Size, const size_t Id>
	__forceinline double Distance2(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		return Length2(u - v);
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Exp(const VecSimd4<Size, Id> & v)
	{
		VecSimd4<Size, Id> result;
		for (int i = 0; i < Size; i++) { result[i] = std::exp(v[i]); }
		return result;
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Floor(const VecSimd4<Size, Id> & v)
	{
		return _mm256_floor_pd(v.v);
	}

	template <const size_t Size, const size_t Id>
	__forceinline bool IsNormalized(const VecSimd4<Size, Id> & a)
	{
		return Math::IsApprox(Length2(a), 1.0);
	}

	template <const size_t Size, const size_t Id>
	__forceinline bool IsZero(const VecSimd4<Size, Id> & a)
	{
		return a == VecSimd4<Size, Id>(0.0);
	}

	template <const size_t Size, const size_t Id>
	__forceinline bool IsFinite(const VecSimd4<Size, Id> & a)
	{
		for (size_t i = 0; i < Size; i++) if (!isfinite2(a[i])) return false;
		return true;
	}

	template <const size_t Size, const size_t Id>
	__forceinline double Length2(const VecSimd4<Size, Id> & v)
	{
		return Dot(v, v);
	}

	template <const size_t Size, const size_t Id>
	__forceinline double Length(const VecSimd4<Size, Id> & v)
	{
		return std::sqrt(Length2(v));
	}

	template <const size_t Size, const size_t Id>
	__forceinline size_t MaxDimension(const VecSimd4<Size, Id> & v)
	{
		size_t index = 0;
		double max = v[0];
		for (size_t i = 1; i < Size; i++) if (v[i] > max) { index = i; max = v[i]; }
		return index;
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Max(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		return _mm256_max_pd(u.v, v.v);
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Min(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		return _mm256_min_pd(u.v, v.v);
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Normalize(const VecSimd4<Size, Id> & v)
	{
		return v / Length(v);
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Pow(const VecSimd4<Size, Id> & v, const double s)
	{
		VecSimd4<Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::pow(v[i], s);
		return result;
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Pow(const VecSimd4<Size, Id> & u, const VecSimd4<Size, Id> & v)
	{
		VecSimd4<Size, Id> result;
		for (size_t i = 0; i < Size; i++) result[i] = std::pow(u[i], v[i]);
		return result;
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Sqrt(const VecSimd4<Size, Id> & v)
	{
		return _mm256_sqrt_pd(v.v);
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Round(const VecSimd4<Size, Id> & v)
	{
        VecSimd4<Size, Id> result;
        for (size_t i = 0; i < Size; i++) result[i] = std::round(v[i]);
        return result;
	}

	template <const size_t Size, const size_t Id>
	__forceinline VecSimd4<Size, Id> Volume(const VecSimd4<Size, Id> & v)
	{
		double result = 1.0;
		for (size_t i = 0; i < Size; i++) result *= v[i];
		return result;
	}
}
