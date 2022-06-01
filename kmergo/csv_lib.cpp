/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2020-07-27
*/

#include "KmerGO.h"

//------------------------------------------------------------
// Read the csv file
// Return true if sucess
//------------------------------------------------------------
bool Read_CSV(std::string csv_path, std::map<std::string, std::string>& stringMap)
{
	std::ifstream input_file(csv_path);
	std::string line;

	if (!input_file)
	{
		return false;
	}

	std::getline(input_file, line);	// check head text
	if (line[line.size() - 1] == '\r')
	{
		line = line.substr(0, line.size() - 1);
	}
	if (line.compare("id,trait") != 0)
	{
		return false;
	}

	while (std::getline(input_file, line))
	{
		if (line[line.size() - 1] == '\r')
		{
			line = line.substr(0, line.size() - 1);
		}
		auto pos = line.find(",");
		stringMap[line.substr(0, pos)] = line.substr(pos + 1);
	}

	return true;
}

