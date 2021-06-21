#pragma once

#include "common/wurst.h"
#include "common/floatimage.h"

namespace Util
{
	struct ErrorImage
	{
		static Fimage<RgbSpectrum> RelDiffRedGreen(const Fimage<RgbSpectrum> & image, const Fimage<RgbSpectrum> & reference)
		{
			assert(image.mSize == reference.mSize);
			Fimage<RgbSpectrum> result(image.mSize);
			result.forEachPixel([&](const Ivec2 & pixel, RgbSpectrum * color)
								{
									RgbSpectrum src1 = image.at(pixel);
									RgbSpectrum src2 = reference.at(pixel);
									const RgbSpectrum d = 2.0 * (src1 - src2) / (Math::Abs(src1 + src2) + RgbSpectrum(1e-3));
									double norm2 = Math::Length2(d);
									double pn = 0.0;
									for (int c = 0; c < RgbSpectrum::NumElements; c++) { pn += d[c]; }
									if (pn > 0.0)
									{
										*color = RgbSpectrum(0.0, norm2, 0.0);
									}
									else if (pn == 0.0)
									{
										*color = RgbSpectrum(0.0, 0.0, 0.0);
									}
									else
									{
										*color = RgbSpectrum(norm2, 0.0, 0.0);
									}
								});
			return result;
		}
	};
}
