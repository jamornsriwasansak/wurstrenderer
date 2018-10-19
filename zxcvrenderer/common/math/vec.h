#pragma once

#include <functional>
#include <cassert>
#include <sstream>

// :( fix ambiguous
inline bool isfinite2(uint64_t) { return true; }
inline bool isfinite2(double x) { return std::isfinite<double>(x); }

template <typename Scalar, const size_t Size, const size_t Id>
struct Vec
{
	friend std::ostream & operator<<(std::ostream & out, const Vec & v)
	{
		out << "(";
		for (size_t i = 0; i < Size - 1; i++)
		{
			out << v[i] << ", ";
		}
		out << v[Size - 1] << ")";
		return out;
	}

	friend Scalar Length2(const Vec & v)
	{
		Scalar sum = 0;
		for (size_t i = 0; i < Size; i++) sum += v[i] * v[i];
		return sum;
	}

	friend Scalar Length(const Vec & v)
	{
		return std::sqrt(Length2(v));
	}

	friend Vec Normalize(const Vec & v)
	{
		return v / Length(v);
	}

	static Vec Abs(const Vec & v)
	{
		Vec result;
		for (size_t i = 0; i < Size; i++) result[i] = std::abs(v[i]);
		return result;
	}

	static Vec Pow(const Vec & v, const Scalar & s)
	{
		Vec result;
		for (size_t i = 0; i < Size; i++) result[i] = std::pow(v[i], s);
		return result;
	}

	static void AssertFinites(const Vec & v)
	{
		for (size_t i = 0; i < Size; i++) assert(isfinite2(v[i]));
	}

	Scalar v[Size];

	Vec()
	{
		for (size_t i = 0; i < Size; i++) v[i] = 0;
		AssertFinites(*this);
	}

	Vec(Scalar s)
	{
		for (size_t i = 0; i < Size; i++) v[i] = s;
		AssertFinites(*this);
	}

	Vec(const Vec & u)
	{
		for (size_t i = 0; i < Size; i++) v[i] = u[i];
		AssertFinites(*this);
	}

	template <typename OtherScalar, size_t Id>
	Vec(const Vec<OtherScalar, Size, Id> & u)
	{
		for (size_t i = 0; i < Size; i++) v[i] = static_cast<Scalar>(u[i]);
		AssertFinites(*this);
	}

	template <typename ... MoreScalars>
	Vec(MoreScalars ... t)
	{
		static_assert(sizeof...(MoreScalars) == Size, "Vec: number of arguments mismatch the type");
		variadicInitialize(0, t...);
		AssertFinites(*this);
	}

	template <typename ... MoreScalars>
	void variadicInitialize(size_t index, Scalar s, MoreScalars ... t)
	{
		v[index] = s;
		variadicInitialize(index + 1, t...);
	}

	void variadicInitialize(size_t index, Scalar s)
	{
		v[index] = s;
	}

	const Scalar & operator[](const size_t index) const
	{
		assert(index < Size);
		return v[index];
	}

	Scalar & operator[](const size_t index)
	{
		assert(index < Size);
		return v[index];
	}

	Vec operate(const Vec & other, const std::function<Scalar(const Scalar a, const Scalar b)>& func) const
	{
		Vec result;
		for (size_t i = 0; i < Size; i++)
		{
			result[i] = func(v[i], other[i]);
		}
		AssertFinites(result);
		return result;
	}

	Vec operator+(const Vec & other) const { return operate(other, [](const Scalar a, const Scalar b) { return a + b; }); }
	Vec operator-(const Vec & other) const { return operate(other, [](const Scalar a, const Scalar b) { return a - b; }); }
	Vec operator*(const Vec & other) const { return operate(other, [](const Scalar a, const Scalar b) { return a * b; }); }
	Vec operator/(const Vec & other) const { return operate(other, [](const Scalar a, const Scalar b) { return a / b; }); }

	Vec operate(const std::function<Scalar(const Scalar a)>& func) const
	{
		Vec result;
		for (size_t i = 0; i < Size; i++)
		{
			result[i] = func(v[i]);
		}
		AssertFinites(result);
		return result;
	}

	Vec operator+(const Scalar other) const { return operate([&](const Scalar a) { return a + other; }); }
	Vec operator-(const Scalar other) const { return operate([&](const Scalar a) { return a - other; }); }
	Vec operator*(const Scalar other) const { return operate([&](const Scalar a) { return a * other; }); }
	Vec operator/(const Scalar other) const { return operate([&](const Scalar a) { return a / other; }); }

	friend Vec operator+(const Scalar other, const Vec & v) { return v + other; }
	friend Vec operator-(const Scalar other, const Vec & v) { return v - other; }
	friend Vec operator*(const Scalar other, const Vec & v) { return v * other; }
	friend Vec operator/(const Scalar other, const Vec & v) { return v.operate([&](const Scalar a) { return other / a; }); }

	Vec & selfOperate(const Vec & other, const std::function<void(Scalar & a, const Scalar & b)>& func)
	{
		for (size_t i = 0; i < Size; i++)
		{
			func(v[i], other[i]);
		}
		AssertFinites(*this);
		return *this;
	}

	Vec & operator+=(const Vec & other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a += b; }); }
	Vec & operator-=(const Vec & other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a -= b; }); }
	Vec & operator*=(const Vec & other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a *= b; }); }
	Vec & operator/=(const Vec & other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a /= b; }); }

	Vec & selfOperate(const Scalar & other, const std::function<void(Scalar & a, const Scalar & b)>& func)
	{
		for (size_t i = 0; i < Size; i++)
		{
			func(v[i], other);
		}
		AssertFinites(*this);
		return *this;
	}

	Vec & operator+=(const Scalar other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a += b; }); }
	Vec & operator-=(const Scalar other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a -= b; }); }
	Vec & operator*=(const Scalar other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a *= b; }); }
	Vec & operator/=(const Scalar other) { return selfOperate(other, [](Scalar & a, const Scalar b) { a /= b; }); }

	bool operator==(const Vec other) const { for (size_t i = 0; i < Size; i++) { if (v[i] != other[i]) return false; } return true; }
	bool operator==(const Scalar other) const { for (size_t i = 0; i < Size; i++) { if (v[i] != other) return false; } return true; }
};
