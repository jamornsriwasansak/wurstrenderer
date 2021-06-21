#pragma once

template <typename Type>
struct KahanSum
{
	KahanSum operator+(const Type & other)
	{
		Type y = other - carry;
		Type t = sum + y;
		KahanSum result;
		result.sum = t;
		result.carry = (t - sum) - y;
		return result;
	}

	KahanSum & operator+=(const Type & other)
	{
		Type y = other - carry;
		Type t = sum + y;
		this->sum = t;
		this->carry = (t - sum) - y;
	}

	Type sum = Type(0.0);
	Type carry = Type(0.0);
};
