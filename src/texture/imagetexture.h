#pragma once

#include "common/wurst.h"

#include "common/floatimage.h"
#include "common/texture.h"

template <typename Type>
struct ImageTexture : public Texture<Type>
{
	ImageTexture(shared_ptr<Fimage<Type>> image): mImage(image)
	{
		mImage->mWrapS = Fimage<Type>::WrapMode::Repeat;
		mImage->mWrapT = Fimage<Type>::WrapMode::Repeat;
	}

	Type eval(const Vec2 & uv) const override
	{
		return mImage->evalBilinear(uv);
	}

	Type average() const override
	{
		return mImage->average();
	}

	shared_ptr<Fimage<Type>> mImage;
};
