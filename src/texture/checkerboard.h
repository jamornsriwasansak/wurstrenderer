#pragma once

#include "common/texture.h"

template <typename Type>
struct CheckerboardTexture : public Texture<Type>
{
	CheckerboardTexture(const Type & color1, const Type & color2, const Uvec2 & numSquares = Uvec2(20, 20)):
		mColor1(color1), mColor2(color2), mNumSquares(numSquares)
	{
	}

	Type eval(const Vec2 & uv) const override
	{
		Uvec2 blockIndex(uv * Vec2(mNumSquares));
		return (blockIndex[0] + blockIndex[1]) % 2 == 0 ? mColor2 : mColor1;
	}

	Type average() const override
	{
		return (mColor1 + mColor2) * 0.5;
	}
	
	Type mColor1;
	Type mColor2;
	Uvec2 mNumSquares;
};