#include "common/wurst.h"
#include "sceneio/simple2parser.h"
#include "sceneio/pbrtparser.h"

void ParsePair(const std::string & option, const std::string & info)
{
	if (option == "--simple2parse")
	{
		Simple2Parser::Parse(filesystem::path(info));
	}
	else if (option == "--jsimple2parse")
	{
		Simple2Parser::Parse(json::parse(info));
	}
}

bool ParseArgs(const int numArgs, const char * args[])
{
	if (numArgs == 1) return false;
	for (int iArg = 1; iArg + 2 <= numArgs; iArg += 2)
	{
		std::string option(args[iArg]);
		std::string info(args[iArg + 1]);
		ParsePair(option, info);
	}
	return true;
}

bool ParseDefaultTxt()
{
	std::ifstream ifs("default.txt");
	if (!ifs.is_open()) { std::ofstream ofs("default.txt"); return false; }

	std::string option, info;
	for (std::string line; std::getline(ifs, line);)
	{
		std::stringstream ss(line);
		ss >> option >> std::quoted(info);
		if (option[0] != '#' && ss.tellg() != -1) ParsePair(option, info);
	}
	return true;
}

int main(const int numArgs, const char * args[])
{
	ParseArgs(numArgs, args) || ParseDefaultTxt();
	return 0;
}
