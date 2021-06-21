#pragma once

#include "common/wurst.h"
#include "ext/threadpool/ThreadPool.h"

struct Parallel
{
	static const size_t DefaultBlock2dSize = 32;

	static int NumThreads()
	{
		// amount of threads in system
		static int numThreads = static_cast<int>(std::thread::hardware_concurrency());
		return numThreads;
	}

	static shared_ptr<ThreadPool> GetThreadPool()
	{
		static shared_ptr<ThreadPool> pool = make_shared<ThreadPool>(NumThreads());
		return pool;
	}

	static int GetThreadIndex()
	{
		auto makeThreadIdMap = [&]()
		{
			std::map<std::thread::id, int> map;
			shared_ptr<ThreadPool> threadPool = GetThreadPool();
			for (Uint i = 0; i < threadPool->workers.size(); i++)
			{
				map[threadPool->workers[i].get_id()] = static_cast<int>(i);
			}
			return map;
		};
		static std::map<std::thread::id, int> threadIdMap = makeThreadIdMap();
		return threadIdMap[std::this_thread::get_id()];
	}

	template <typename Scalar, typename Func>
	static void Split(const Scalar start, const Scalar end, const Scalar chunkSize, const Func & func)
	{
		shared_ptr<ThreadPool> pool = GetThreadPool();
		std::vector<std::future<void>> futures;
		for (Scalar i = start; i < end; i += chunkSize)
		{
			Scalar taskStart = i;
			Scalar taskEnd = std::min(i + chunkSize, end);
			futures.emplace_back(pool->enqueue(func, taskStart, taskEnd));
		}
		for (std::future<void> & future : futures) { future.get(); }
	}

	template <typename Scalar, typename Func>
	static void Split(const Scalar start, const Scalar end, const Func & func)
	{
		Split(start, end, static_cast<Scalar>((end - start) / (2 * NumThreads()) + 1), func);
	}

	template <typename Scalar, typename Func>
	static void ForEach(const Scalar start, const Scalar end, const Func & func)
	{
		Split(start, end, [&](const Scalar a, const Scalar b) { for (Scalar i = a; i < b; i++) { func(i); } });
	}

	static inline void Split2d(const Ivec2 & start, const Ivec2 & end, const int maxBlockSize, const std::function<void(const Ibound2 &)> & func)
	{
		std::vector<std::thread> workers;
		std::mutex mutex;

		std::vector<std::pair<const Ivec2, const Ivec2>> blocks;

		// generate image blocks
		std::function<void(const Ivec2 & start, const Ivec2 & end)> generateBlocks;
		generateBlocks = [&](const Ivec2 & start, const Ivec2 & end)
		{
			const Ivec2 blockSize = end - start;
			const int maxExtent = std::max(blockSize[0], blockSize[1]);
			if (maxExtent <= maxBlockSize) // default block 2d size
			{
				const int minExtent = std::min(blockSize[0], blockSize[1]);
				if (minExtent != 0)
				{
					// check if block is not degenerate
					blocks.emplace_back(std::pair<const Ivec2, const Ivec2>(start, end));
				}
			}
			else if (maxExtent == blockSize[0])
			{
				// split along x axis
				const int halfNumBlocksX = Math::DivAndCeil(blockSize[0], maxBlockSize) / 2;
				const int halfX = start[0] + std::min(halfNumBlocksX * maxBlockSize, blockSize[0]);
				generateBlocks(start, Ivec2(halfX, end[1]));
				generateBlocks(Ivec2(halfX, start[1]), end);
			}
			else if (maxExtent == blockSize[1])
			{
				// split along y axis
				const int halfNumBlocksY = Math::DivAndCeil(blockSize[1], maxBlockSize) / 2;
				const int halfY = start[1] + std::min(halfNumBlocksY * maxBlockSize, blockSize[1]);
				generateBlocks(start, Ivec2(end[0], halfY));
				generateBlocks(Ivec2(start[0], halfY), end);
			}
		};
		generateBlocks(start, end);

		shared_ptr<ThreadPool> pool = GetThreadPool();
		std::vector<std::future<void>> futures;
		for (std::pair<const Ivec2, const Ivec2> block : blocks)
		{
			Ibound2 b(block.first, block.second);
			futures.emplace_back(pool->enqueue(func, b));
		}

		for (std::future<void> & future : futures) { future.get(); }
	}

	static inline void Split2d(const Ivec2 & size, const std::function<void(const Ibound2 &)> & func)
	{
		Split2d(Ivec2(0, 0), Ivec2(size), DefaultBlock2dSize, func);
	}

	static inline void ForEach2d(const Ivec2 & size, const std::function<void(const Ivec2 &)> & func)
	{
		Split2d(Ivec2(0, 0), Ivec2(size), DefaultBlock2dSize, [&](const Ibound2 & block)
		{
			for (const Ivec2 p : block) { func(p); }
		});
	}
};


struct Serial
{
	static inline void Split2d(const Ivec2 & start, const Ivec2 & end, const int maxBlockSize, const std::function<void(const Ibound2 &)> & func)
	{
		std::vector<std::pair<const Ivec2, const Ivec2>> blocks;

		// generate image blocks
		std::function<void(const Ivec2 & start, const Ivec2 & end)> generateBlocks;
		generateBlocks = [&](const Ivec2 & start, const Ivec2 & end)
		{
			const Ivec2 blockSize = end - start;
			const int maxExtent = std::max(blockSize[0], blockSize[1]);
			if (maxExtent <= maxBlockSize) // 32 pixels
			{
				blocks.emplace_back(std::pair<const Ivec2, const Ivec2>(start, end));
			}
			else if (maxExtent == blockSize[0])
			{
				// split along x axis
				const int halfX = (start[0] + end[0]) / 2;
				generateBlocks(start, Ivec2(halfX, end[1]));
				generateBlocks(Ivec2(halfX, start[1]), end);
			}
			else if (maxExtent == blockSize[1])
			{
				// split along y axis
				const int halfY = (start[1] + end[1]) / 2;
				generateBlocks(start, Ivec2(end[0], halfY));
				generateBlocks(Ivec2(start[0], halfY), end);
			}
		};
		generateBlocks(start, end);

		for (std::pair<const Ivec2, const Ivec2> block : blocks)
		{
			func(Ibound2(block.first, block.second));
		}
	}

	template <typename Scalar, typename Func>
	static void Split(const Scalar start, const Scalar end, const Scalar chunkSize, const Func & func)
	{
		for (Scalar i = start; i < end; i += chunkSize) { func(start, end); }
	}

	template <typename Scalar, typename Func>
	static void Split(const Scalar start, const Scalar end, const Func & func)
	{
		Split(start, end, ((end - start) / (2 * Parallel::NumThreads()) + 1), func);
	}

	template <typename Scalar, typename Func>
	static void ForEach(const Scalar start, const Scalar end, const Func & func)
	{
		Split(start, end, [&](const Scalar a, const Scalar b) { for (Scalar i = a; i < b; i++) { func(i); } });
	}

	static inline void Split2d(const Ivec2 & size, const std::function<void(const Ibound2 &)> & func)
	{
		Split2d(Ivec2(0, 0), size, Parallel::DefaultBlock2dSize, func);
	}

	static inline void ForEach2d(const Ivec2 & size, const std::function<void(const Ivec2 &)> & func)
	{
		Split2d(Ivec2(0, 0), size, Parallel::DefaultBlock2dSize, [&](const Ibound2 & bound) {
			for (int y = bound.pMin[1]; y < bound.pMax[1]; y++)
				for (int x = bound.pMin[0]; x < bound.pMax[0]; x++)
					func(Ivec2(x, y));
		});
	}
};