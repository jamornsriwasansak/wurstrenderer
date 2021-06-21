#pragma once

#include <common/wurst.h>

template <typename Object>
struct HashGrid
{
	int hash3(const Ivec3 &p, int hashSize) const
	{
		return (p[0] * mResolution[0] * mResolution[1] + p[1] * mResolution[1] + p[2]) % hashSize;
	}

	HashGrid(const Bound3 & sceneBound, const Ivec3 & resolution, const int hashSize):
		mBound(sceneBound),
		mResolution(resolution),
		mHash(hashSize),
		mHashMutex(hashSize)
	{
	}

	void insert(const Vec3 & position, const double radius, const Object & object)
	{
		Vec3 radius3(radius);
		Ivec3 grid0 = toGrid(position - radius3);
		Ivec3 grid1 = toGrid(position + radius3);

		for (int z = grid0[2]; z <= grid1[2]; z++)
			for (int y = grid0[1]; y <= grid1[1]; y++)
				for (int x = grid0[0]; x <= grid1[0]; x++)
				{
					const int hashValue = hash3(Ivec3(x, y, z), int(mHash.size()));
					// acquire lock
					{
						std::lock_guard<std::mutex> lock(mHashMutex[hashValue]);
						mHash[hashValue].emplace_back(object);
					}
				}
	}

	const std::vector<Object> & getSearchBucket(const Vec3 & position) const
	{
		// note no lock during get. leads to bug when the program inserts into and gets bucket at the same time
		const int hashValue = hash3(toGrid(position), int(mHash.size()));
		return mHash[hashValue];
	}

	Ivec3 toGrid(const Vec3 & position) const 
	{
		const Vec3 ratioPos = (position - mBound.pMin) / (mBound.pMax - mBound.pMin);
		assert(Bound3::IsInside(Bound3(Vec3(0.0), Vec3(1.0)), ratioPos));
		const Ivec3 gridPos = Ivec3(ratioPos * Vec3(mResolution));
		return gridPos;
	}

	Bound3 mBound;
	Ivec3 mResolution;
	std::vector<std::mutex> mHashMutex;
	std::vector<std::vector<Object>> mHash;
};
