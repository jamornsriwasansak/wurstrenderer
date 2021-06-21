#pragma once

#include "common/wurst.h"

namespace Util
{
	struct Json
	{
		static Vec3 Vec3FromJson(const json & parentJson, const std::string & key)
		{
			const json json = parentJson[key];
			if (json.size() == 1) return Vec3(static_cast<double>(json));
			if (json.size() == 3) return Vec3(json[0], json[1], json[2]);
			throw std::runtime_error((std::string("can't parse ") + key).c_str());
		}

		static Vec2 Vec2FromJson(const json & parentJson, const std::string & key)
		{
			const json json = parentJson[key];
			if (json.size() == 1) return Vec2(static_cast<double>(json));
			if (json.size() == 2) return Vec2(json[0], json[1]);
			throw std::runtime_error((std::string("can't parse ") + key).c_str());
		}

		static RgbSpectrum RgbFromJson(const json & parentJson, const std::string & key)
		{
			return SpectrumConvertUtil::SpectrumFromRgb(Vec3FromJson(parentJson, key));
		}

		static Uvec2 Uvec2FromJson(const json & parentJson, const std::string & key)
		{
			const json json = parentJson[key];
			if (json.size() != 2) throw std::runtime_error((std::string("can't parse ") + key).c_str());
			return Uvec2(json[0], json[1]);
		}

		static Ivec2 Ivec2FromJson(const json & parentJson, const std::string & key)
		{
			const json json = parentJson[key];
			if (json.size() != 2) throw std::runtime_error((std::string("can't parse ") + key).c_str());
			return Ivec2(json[0], json[1]);
		}

		template <typename T>
		static T GetValue(const json & json, const std::string & key, const T defaultValue)
		{
			if (json.find(key) == json.end()) { return defaultValue; }
			return json[key];
		}
	};
}