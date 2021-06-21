#pragma once

#include <common/wurst.h>

//template <typename Object>

using Object = int;

struct BvhNnNode
{
	BvhNnNode()
	{
	}

	Bound3 mBound;
};

struct BvhNnSearch
{

	BvhNnSearch()
	{
	}

	void construct(std::vector<std::pair<Vec3, Object>> & objects, const int start, const int end)
	{
		const int numObjects = end - start;

		// constructing bound of all positions
		Bound3 bound;
		for (const std::pair<Vec3, Object> & object : objects) { bound = Bound3::Union(bound, object.first); }

		// find the longest axis
		const Vec3 boundSize = bound.pMax - bound.pMin;
		const size_t maxExtentDim = Math::MaxDimension(boundSize);

		// sort the objects along the longest axis
		std::sort(objects.begin() + start, objects.begin() + end,
				  [&](const std::pair<Vec3, Object> & object1, const std::pair<Vec3, Object> & object2) { return object1.first[maxExtentDim] < object2.first[maxExtentDim]; });

		// compute surface area from left and from right
		Bound3 bLeft;
		std::vector<double> surfaceAreaFromLeft(end - start);
		for (int iObject = start; iObject < end; iObject++)
		{
			bLeft = Bound3::Union(bLeft, objects[iObject].first);
			surfaceAreaFromLeft[iObject - start] = Bound3::Area(bLeft);
		}

		Bound3 bRight;
		std::vector<double> surfaceAreaFromRight(end - start);
		for (int iObject = end - 1; iObject >= start; iObject--)
		{
			bRight = Bound3::Union(bRight, objects[iObject].first);
			surfaceAreaFromRight[iObject - start] = Bound3::Area(bRight);
		}

		// find index that yield minimum surface area
		double minSurfaceArea = std::numeric_limits<double>::infinity();
		int minSurfaceAreaIndex = 0;
		for (int i = 0; i < numObjects; i++)
		{
			const double surfaceArea = double(i) * surfaceAreaFromLeft[i] + double(numObjects - i - 1) * surfaceAreaFromRight[i + 1];
			if (surfaceArea < minSurfaceArea)
			{
				minSurfaceAreaIndex = i;
				minSurfaceArea = surfaceArea;
			}
		}

		// recursively split
		construct(objects, start, minSurfaceAreaIndex + 1);
		construct(objects, minSurfaceAreaIndex + 1, end);
	}
};
