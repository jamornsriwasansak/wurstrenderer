#pragma once

#include "common/wurst.h"
#include "common/util/file.h"

namespace Util
{
	struct Vector
	{
		template <typename Type>
		static std::vector<std::vector<Type>> CreateVector2d(const int numRows, const int numCols)
		{
			std::vector<std::vector<Type>> result(numRows);
			for (int i = 0; i < numRows; i++) result[i] = std::vector<Type>(numCols);
			return result;
		}

		template <typename Type>
		static std::vector<std::vector<Type>> CreateVector2d(const std::vector<int> & lengths)
		{
			std::vector<std::vector<Type>> result(lengths.size());
			for (int i = 0; i < lengths.size(); i++) result[i] = std::vector<Type>(lengths[i]);
			return result;
		}

		template <typename Type>
		static std::vector<Type> LoadVectorBin(const filesystem::path & filepath)
		{
			auto filesize = filesystem::file_size(filepath);
			std::ifstream is(filepath, std::ios::in | std::ios::binary);
			std::vector<Type> vec(filesize / sizeof(Type));
			is.read((char*)(&vec[0]), filesize);

			return vec;
		}

		template <typename Type>
		static void WriteVectorBin(const filesystem::path & filepath, const std::vector<Type> & datas)
		{
			Util::File::SafeWritePath(filepath);
			std::ofstream os(filepath, std::ios::out | std::ios::binary);
			os.write((char*)&datas[0], datas.size() * sizeof(Type));
			os.close();
		}

		template <typename Type>
		static void WriteVector2dBin(const filesystem::path & filepath, const std::vector<std::vector<Type>> & datas)
		{
			Util::File::SafeWritePath(filepath);
			std::ofstream os(filepath, std::ios::binary | std::ios::out);
			std::vector<uint64_t> sizes;
			for (int i = 0; i < datas.size(); i++)
			{
				throw std::runtime_error("unimpl");
			}
		}
	};
};
