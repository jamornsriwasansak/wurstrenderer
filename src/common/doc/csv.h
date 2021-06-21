#pragma once

#include "common/wurst.h"

struct CsvUtil
{
	template <typename Type>
	static void Save(std::vector<Type> & datas, const filesystem::path & filepath)
	{
		std::ofstream of(filepath);
		if (!of.is_open()) throw std::runtime_error("can't open filepath");
		for (Uint i = 0; i < datas.size(); i++)
		{
			of << datas[i];
			if (i + 1 != datas.size()) of << ", ";
		}
		of.close();
	}
};