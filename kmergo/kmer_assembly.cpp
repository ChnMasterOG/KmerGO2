/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2022-05-30
*/

#include "KmerGO.h"
#include "KmerGO_io.h"

//--------------------------------------------------------------------------------------------------------------------
short kmer_assembly(KmerGO_parameters& p, std::string& loginfo)
{
	std::string onput_file_name_A = joinPath(p.results_path, "Afeatures.fa");
	std::string onput_file_name_B = joinPath(p.results_path, "Bfeatures.fa");
	std::ofstream output_file_A(onput_file_name_A);
	std::ofstream output_file_B(onput_file_name_B);
	std::vector<std::string> files_vec;
	std::string line;
	std::string NameA, NameB;
	uint64 countA = 0, countB = 0;

	getFiles(p.features_path, files_vec, false);

	for (auto f : files_vec)
	{
		if (f.compare(0, 11, "categorical") == 0)
		{
			if (p.mode != 0)
			{
				loginfo = "Mode Error!";
				return 1;
			}
		}
		else if (f.compare(0, 10, "continuous") == 0)
		{
			if (p.mode != 1)
			{
				loginfo = "Mode Error!";
				return 1;
			}
		}
		else
		{
			loginfo = "Unrecognized Files!";
			return 1;
		}
		std::ifstream input_file(joinPath(p.features_path, f));
		std::getline(input_file, line);	// skip headtext
		while (std::getline(input_file, line))
		{
			if (line[line.size() - 1] == '\r')
			{
				line = line.substr(0, line.size() - 1);
			}
			std::string label;
			if (p.mode == 0)
			{
				label = line.substr(line.rfind('\t') + 1);
			}
			else
			{
				label = "Result";
			}
			if (NameA.empty())
			{
				NameA = label;
			}
			else if (NameB.empty() && NameA.compare(label) != 0)
			{
				NameB = label;
			}
			if (NameA.compare(label) == 0)
			{
				countA++;
				output_file_A << '>' << std::to_string(countA) << '\n' << line.substr(0, line.find('\t')) << '\n';
			}
			else if (NameB.compare(label) == 0)
			{
				countB++;
				output_file_B << '>' << std::to_string(countB) << '\n' << line.substr(0, line.find('\t')) << '\n';
			}
		}
	}

	uint32 r, r1 = 0, r2 = 0;
	std::string cmd, current_path;
	current_path = getCWD();

	if (!NameA.empty())
	{
		cmd = "mv " + onput_file_name_A + " " + joinPath(p.results_path, NameA + "_specific.fa");
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "mv Command Error!";
			return 3;
		}
	}
	else
	{
		cmd = "rm " + joinPath(p.results_path, "Afeatures.fa");
		r = system(cmd.c_str());
	}
	if (!NameB.empty())
	{
		cmd = "mv " + onput_file_name_B + " " + joinPath(p.results_path, NameB + "_specific.fa");
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "mv Command Error!";
			return 3;
		}
	}
	else
	{
		cmd = "rm " + joinPath(p.results_path, "Bfeatures.fa");
		r = system(cmd.c_str());
	}

	if (!NameA.empty())
	{
		cmd = joinPath(current_path, "bin/cap3 ") + joinPath(p.results_path, NameA + "_specific.fa") + " -i 30  -j 31  -o 18  -s 300";
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "CAP3 Software Error!";
			return 2;
		}
		cmd = "mv " + joinPath(p.results_path, NameA + "_specific.fa") + " " + joinPath(p.results_path, NameA + "_specific_kmer.fa");
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "mv Command Error!";
			return 3;
		}
	}
	else
	{
		r1 = 1;
	}
	
	if (!NameB.empty())
	{
		cmd = joinPath(current_path, "bin/cap3 ") + joinPath(p.results_path, NameB + "_specific.fa") + " -i 30  -j 31  -o 18  -s 300";
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "CAP3 Software Error!";
			return 2;
		}
		cmd = "mv " + joinPath(p.results_path, NameB + "_specific.fa") + " " + joinPath(p.results_path, NameB + "_specific_kmer.fa");
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "mv Command Error!";
			return 3;
		}
	}
	else
	{
		r2 = 1;
	}
	
	if (r1 && r2)
	{
		loginfo = "No specific k-mer!";
		return 5;
	}

	return 0;
}


