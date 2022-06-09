/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2022-05-30
*/

#include "KmerGO.h"
#include "KmerGO_io.h"
#include "../kmc_api/kmc_file.h"

KmerGO_tools tool = NO_SELECTED;
KmerGO_parameters param;
std::string loginfo;

// -------------------------------------------------------------------------
// Check if --help or --version was used
// -------------------------------------------------------------------------
bool help_or_version(int argc, char** argv)
{
	const std::string version = "--version";
	const std::string help = "--help";
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == version || argv[i] == help)
			return true;
	}
	return false;
}

// -------------------------------------------------------------------------
// Print execution options 
// -------------------------------------------------------------------------
void print_info(void)
{
	std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
		<< "\nUsage:\nKmerGO <tool> [options]\n"
		<< "Tool:\n"
		<< "one of kmc3, union, filtering, cap3 or all\n"
		<< "Options:\n"
		<< "-i <string>: input files path (required)\n"
		<< "-o <string>: output files path (required)\n"
		<< "-t <string>: a csv file path of trait information (required when tool is filtering or all)\n"
		<< "-e <string>: temp files path (required when tool is kmc3 or all, default: ./)\n"
		<< "-m <int>: mode, 0-categorical, 1-continuous (default: 0)\n"
		<< "-k <int>: k-mer length (k from 14 to 256, default: 40)\n"
		<< "-ci <int>: minimal k-mer occurring times (default: 2)\n"
		<< "-cs <int>: maximal k-mer occurring times (default: 65535)\n"
		<< "-n <int>: number of processes (n from 1 to 256, default: 24)\n"
		<< "-assl <float>: when mode = 0, logical features ASS value (default: 0.8)\n"
		<< "-p <float>: numeric(mode=0) or logical(mode=1) features rank sum test p threshold value (default: 0.01)\n"
		<< "-assn <float>: when mode = 0, numeric features logistic regression ASS value (default: 0.8)\n"
		<< "-corr <float>: when mode = 1, numeric features coefficient of association threshold value (default: 0.8)\n";
}

// -------------------------------------------------------------------------
// Main function 
// -------------------------------------------------------------------------
int main(int argc, char* argv[])
{

	if (argc <= 2 || help_or_version(argc, argv))
	{
		print_info();
		return EXIT_FAILURE;
	}

	int32 i;
	char char_temp[1024];
	std::string input_path = "No select", output_path = "No select";

	//------------------------------------------------------------
	// Parse input parameters
	//------------------------------------------------------------

	if (strncmp(argv[1], "kmc3", 4) == 0)
		tool = KMC3;
	else if (strncmp(argv[1], "union", 5) == 0)
		tool = UNION;
	else if (strncmp(argv[1], "filtering", 9) == 0)
		tool = FILTERING;
	else if (strncmp(argv[1], "cap3", 4) == 0)
		tool = CAP3;
	else if (strncmp(argv[1], "all", 3) == 0)
	{
		param.single_step = false;
		tool = ALL;
	}

	if (tool == NO_SELECTED)
	{
		std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
			<< "Error! <tool> is required and is one of {kmc3,union,filtering,cap3,all}.\n";
		return EXIT_FAILURE;
	}

	for (i = 2; i < argc; ++i)
	{
		if (argv[i][0] == '-')
		{
			if (strncmp(argv[i], "-m", 2) == 0)
				param.mode = atoi(&argv[i][3]);
			else if (strncmp(argv[i], "-k", 2) == 0)
				param.k_value = atoi(&argv[i][3]);
			else if (strncmp(argv[i], "-ci", 3) == 0)
				param.ci_value = atoi(&argv[i][4]);
			else if (strncmp(argv[i], "-cs", 3) == 0)
				param.cs_value = atoi(&argv[i][4]);
			else if (strncmp(argv[i], "-n", 2) == 0)
				param.process_number = atoi(&argv[i][3]);
			else if (strncmp(argv[i], "-p", 2) == 0)
				param.p_value = atof(&argv[i][3]);
			else if (strncmp(argv[i], "-assl", 5) == 0)
				param.assl = atof(&argv[i][6]);
			else if (strncmp(argv[i], "-assn", 5) == 0)
				param.assn = atof(&argv[i][6]);
			else if (strncmp(argv[i], "-corr", 5) == 0)
				param.corr_vaule = atof(&argv[i][6]);
			else if (strncmp(argv[i], "-i", 2) == 0)
			{
				strcpy(char_temp, &argv[i][3]);
				input_path = char_temp;
				//if (input_path.compare(input_path.length() - 1, 1, "/", 0, 1) != 0) input_path = input_path + "/";
			}
			else if (strncmp(argv[i], "-o", 2) == 0)
			{
				strcpy(char_temp, &argv[i][3]);
				output_path = char_temp;
				//if (output_path.compare(output_path.length() - 1, 1, "/", 0, 1) != 0) output_path = output_path + "/";
			}
			else if (strncmp(argv[i], "-t", 2) == 0)
			{
				strcpy(char_temp, &argv[i][3]);
				param.csv_path = char_temp;
			}
			else if (strncmp(argv[i], "-e", 2) == 0)
			{
				strcpy(char_temp, &argv[i][3]);
				param.temp_path = char_temp;
				//if (param.temp_path.compare(param.temp_path.length() - 1, 1, "/", 0, 1) != 0) param.temp_path = param.temp_path + "/";
			}
			else
			{
				std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
					<< "Error! Can not read some parameters.\n";
				return EXIT_FAILURE;
			}
		}
	}

	if (input_path == "No select" or output_path == "No select")
	{
		std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
			<< "Error! input files path and output files path are required.\n";
		return EXIT_FAILURE;
	}

	switch (tool)
	{
	case KMC3:
		param.samples_path = input_path;
		param.kmercountings_path = output_path;
		break;
	case UNION:
		param.kmercountings_path = input_path;
		param.matrixes_path = output_path;
		break;
	case FILTERING:
		param.matrixes_path = input_path;
		param.features_path = output_path;
		break;
	case CAP3:
		param.features_path = input_path;
		param.results_path = output_path;
		break;
	case ALL:
		param.samples_path = input_path;
		param.kmercountings_path = joinPath(param.temp_path, "kmer_countings/");
		param.matrixes_path = joinPath(param.temp_path, "kmer_matrixes/");
		param.features_path = joinPath(param.temp_path, "kmer_features/");
		param.results_path = output_path;
		// create folder
		if (system(("mkdir " + param.kmercountings_path).c_str()))
		{
			std::cout << "Error: mkdir" << std::endl;
		}
		if (system(("mkdir " + param.matrixes_path).c_str()))
		{
			std::cout << "Error: mkdir" << std::endl;
		}
		if (system(("mkdir " + param.features_path).c_str()))
		{
			std::cout << "Error: mkdir" << std::endl;
		}
	case NO_SELECTED:
		;
	}

	//------------------------------------------------------------
	// Read CSV
	//------------------------------------------------------------

	bool csv_err = false;
	if (tool == FILTERING || tool == ALL)
	{
		std::cout << "read CSV..." << "\n";
		csv_err = !Read_CSV(param.csv_path, param.trait_information_map);
		if (csv_err)
		{
			std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
				<< "Error! Can't not read the csv file.\n";
			return EXIT_FAILURE;
		}
		if (param.trait_information_map.size() <= 1)
		{
			std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
				<< "Error! Sample number is too small.\n";
			return EXIT_FAILURE;
		}
		if (param.mode == 1)
		{
			param.groupA_name = "";
			param.groupA_number = param.trait_information_map.size();
		}
		else
		{
			param.groupA_name = param.trait_information_map.begin()->second;
			param.groupB_name = "";
			std::map<std::string, std::string>::iterator it;
			for (it = param.trait_information_map.begin(); it != param.trait_information_map.end(); it++) {
				if (!param.groupA_name.compare(it->second)) param.groupA_number++;
				else if (param.groupB_name.compare(it->second) && param.groupB_name.empty())
				{
					param.groupB_name = it->second;
					param.groupB_number++;
				}
				else if (!param.groupB_name.compare(it->second)) param.groupB_number++;
				else 
				{
					std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
						<< "Error! More than 2 groups of data.\n";
					return EXIT_FAILURE;
				}
			}
		}
	}

	//------------------------------------------------------------
	// Run KmerGO
	//------------------------------------------------------------

	///////////////////////////////////////////////////////
	////////////////////////Module 1///////////////////////
	//////////////////////////KMC3/////////////////////////
	///////////////////////////////////////////////////////
	if (tool == KMC3 || tool == ALL)
	{
		std::cout << "runing KMC3..." << "\n";
		if (kmer_counting(param, loginfo))
		{
			std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
				<< "KMC3 step error! \n" << loginfo << "\n";
			return EXIT_FAILURE;
		}
	}

	///////////////////////////////////////////////////////
	////////////////////////Module 2///////////////////////
	/////////////////////////Union/////////////////////////
	///////////////////////////////////////////////////////
	if (tool == UNION || tool == ALL)
	{
		std::cout << "runing UNION..." << "\n";
		if (kmer_union(param, loginfo))
		{
			std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
				<< "Union step error! \n" << loginfo << "\n";
			return EXIT_FAILURE;
		}
	}

	////////////////////////Module 5///////////////////////
	/////////////////Add PCA function here/////////////////
	/////////Read k-mer matrix and calculate PCA///////////
	/* Add code here. */

	///////////////////////////////////////////////////////
	////////////////////////Module 3///////////////////////
	///////////////////////Filtering///////////////////////
	///////////////////////////////////////////////////////
	if (tool == FILTERING || tool == ALL)
	{
		std::cout << "runing FILTERING..." << "\n";
		if (kmer_filtering(param, loginfo))
		{
			std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
				<< "Filtering step error! \n" << loginfo << "\n";
			return EXIT_FAILURE;
		}
	}

	///////////////////////////////////////////////////////
	////////////////////////Module 4///////////////////////
	//////////////////////////CAP3/////////////////////////
	///////////////////////////////////////////////////////
	if (tool == CAP3 || tool == ALL)
	{
		std::cout << "runing CAP3..." << "\n";
		if (kmer_assembly(param, loginfo))
		{
			std::cout << "KmerGO ver. " << KMERGO_VER << " (" << KMERGO_DATE << ")\n"
				<< "CAP3 step error! \n" << loginfo << "\n";
			return EXIT_FAILURE;
		}
	}

	return EXIT_SUCCESS;
}

