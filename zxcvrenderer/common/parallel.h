#pragma once

#include "common/asdfjkl.h"
#include <iostream>
#include <memory>
#include <vector>
#include <mutex>
#include <thread>
#include <atomic>
#include <algorithm>

struct Parallel
{
	static inline uint32_t NumThreads() // amount of threads in system
	{
		static uint32_t numThreads = std::thread::hardware_concurrency();
		return numThreads;
	}

	static void Split(const int start, const int end, const int chunkSizePerThread, const std::function<void(int, int)> & func)
	{
		static uint32_t numThreads = NumThreads();
		std::vector<std::thread> workers;
		std::mutex mutex;

		int iter = start;
		auto iterate = [&](int * newStart, int * newEnd)
		{
			std::unique_lock<std::mutex> lock(mutex);
			if (iter >= end) { return false; }
			*newStart = iter;
			*newEnd = std::min(iter + chunkSizePerThread, end);
			iter += chunkSizePerThread;
			return true;
		};

		for (uint32_t i = 0; i < numThreads; i++)
		{
			workers.push_back(std::thread([&]() {
				int newStart, newEnd;
				while (iterate(&newStart, &newEnd)) { func(newStart, newEnd); }
			}));
		}
		for (std::thread & worker : workers) { worker.join(); }
	}

	static void Split(const int start, const int end, const std::function<void(int, int)> & func)
	{
		static uint32_t numThreads = NumThreads();
		const int chunkSizePerThread = (end - start + numThreads - 1) / numThreads;
		Split(start, end, chunkSizePerThread, func);
	}

	static void Split2d(const Uvec2 & start, const Uvec2 & end, const uint32_t maxBlockSize, const std::function<void(const Uvec2 &, const Uvec2 &)> & func)
	{
		static uint32_t numThreads = NumThreads();
		std::vector<std::thread> workers;
		std::mutex mutex;

		std::vector<std::pair<const Uvec2, const Uvec2>> blocks;

		// generate image blocks
		std::function<void(const Uvec2 & start, const Uvec2 & end)> generateBlocks;
		generateBlocks = [&](const Uvec2 & start, const Uvec2 & end)
		{
			const Uvec2 blockSize = end - start;
			const uint32_t maxExtent = std::max(blockSize[0], blockSize[1]);
			if (maxExtent <= maxBlockSize) // 32 pixels
			{
				blocks.emplace_back(std::pair<const Uvec2, const Uvec2>(start, end));
			}
			else if (maxExtent == blockSize[0])
			{
				// split along x axis
				const uint32_t halfX = (start[0] + end[0]) / 2;
				generateBlocks(start, Uvec2(halfX, end[1]));
				generateBlocks(Uvec2(halfX, start[1]), end);
			}
			else if (maxExtent == blockSize[1])
			{
				// split along y axis
				const size_t halfY = (start[1] + end[1]) / 2;
				generateBlocks(start, Uvec2(end[0], halfY));
				generateBlocks(Uvec2(start[0], halfY), end);
			}
		};

		generateBlocks(start, end);
		size_t iter = 0;

		auto iterate = [&](Uvec2 * newStart, Uvec2 * newEnd)
		{
			std::unique_lock<std::mutex> lock(mutex);
			if (iter >= blocks.size()) { return false; }
			*newStart = blocks[iter].first;
			*newEnd = blocks[iter].second;
			iter++;
			return true;
		};

		for (uint32_t i = 0; i < numThreads; i++)
		{
			workers.push_back(std::thread([&]() { Uvec2 newStart, newEnd; while (iterate(&newStart, &newEnd)) { func(newStart, newEnd); } }));
		}
		for (std::thread & worker : workers) { worker.join(); }
	}
};
