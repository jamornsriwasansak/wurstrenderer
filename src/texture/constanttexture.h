#pragma once

#include "common/texture.h"

template <typename Type>
struct ConstantTexture : public Texture<Type>
{
	ConstantTexture(const Type & value): mValue(value) {}

	Type eval(const Vec2 &) const override { return mValue; }

	Type average() const override { return mValue; }

	Type mValue;
};
