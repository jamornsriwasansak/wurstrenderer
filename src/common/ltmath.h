#pragma once

#include "common/math/math.h"

// light transport mathematics
struct LtMath
{
    static double ApproximateFresnelTerm(const Vec3 & i, const Vec3 & n, const double relativeRefractionIndex)
    {
        // optical properties of inhomogeneous meterials, egan et al.
        return -1.440 / (relativeRefractionIndex * relativeRefractionIndex) + 0.710 / relativeRefractionIndex + 0.668 + 0.0636 * relativeRefractionIndex;
    }
};
