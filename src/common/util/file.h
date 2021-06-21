#pragma once

#include "common/wurst.h"

namespace Util
{
	struct File
	{
		static Uint CountFiles(const filesystem::path & filepath, const std::string & extension = "!")
		{
			Uint result = 0;
			for (filesystem::directory_iterator it(filepath); it != filesystem::directory_iterator(); ++it)
			{
				if (filesystem::is_directory(it->path())) { continue; }
				if (it->path().extension().string() == extension || extension == "!") { result += 1; }
			}
			return result;
		}

		static std::vector<filesystem::path> GetFilepathsInFolder(const filesystem::path & filepath, const std::vector<std::string> & extensions = { "!" })
		{
			std::vector<filesystem::path> paths;
			for (filesystem::directory_iterator it(filepath); it != filesystem::directory_iterator(); ++it)
			{
				if (filesystem::is_directory(it->path())) { continue; }
				for (const std::string & extension : extensions)
					if (it->path().extension().string() == extension || extension == "!") { paths.emplace_back(it->path()); break; }
			}
			return paths;
		}

		static json JsonArrayFromVec3(const Vec3 & v)
		{
			json result = json::array();
			for (Uint i = 0; i < 3; i++) { result.push_back(v[i]); }
			return result;
		}

		static filesystem::path MakeAbsolutePath(const filesystem::path & path)
		{
			if (path.is_absolute()) return path;
			return filesystem::current_path() / path;
		}

		static filesystem::path SafeOverwritePath(const filesystem::path & filepath)
		{
			// create folder if not exist
			filesystem::create_directories(filepath.parent_path());

			filesystem::path targetFilepath = filepath;
			filesystem::path parentPath = filepath.parent_path();
			filesystem::path filename = filepath.stem().filename();
			filesystem::path extension = filepath.extension();
			return targetFilepath;
		}

		static filesystem::path SafeWritePath(const filesystem::path & filepath)
		{
			// create folder if not exist
			filesystem::create_directories(filepath.parent_path());

			filesystem::path targetFilepath = filepath;
			filesystem::path parentPath = filepath.parent_path();
			filesystem::path filename = filepath.stem().filename();
			filesystem::path extension = filepath.extension();
			Uint counter = 0;
			while (filesystem::exists(targetFilepath))
			{
				counter += 1;
				const std::string mark = filename.string() + "(" + std::to_string(counter) + ")" + extension.string();
				targetFilepath = parentPath / mark;
			}
			return targetFilepath;
		}

		static void SafeWriteString(const filesystem::path & filepath, const std::string & str)
		{
			SafeWritePath(filepath);
			std::ofstream os(filepath);
			if (!os.is_open()) throw std::runtime_error("cannot open file");
			os << str;
			os.close();
		}
	};
}