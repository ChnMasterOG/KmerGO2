/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2022-05-30
*/

#include "KmerGO.h"
#include "KmerGO_io.h"

static volatile uint64 progress = 0;
static volatile short filtering_err_number = 0;
static std::mutex mtx;

//--------------------------------------------------------------------------------------------------------------------
static void show_Progress(uint64 allfilesize)
{
	printf("STEP k-mer filtering start!\n");
	while (progress < allfilesize)
	{
		printf("\rTotal: %.2f%%  ", (double)progress * 100 / allfilesize);
	}
	printf("\rTotal: 100.00%%\n");
}

//--------------------------------------------------------------------------------------------------------------------
void categorical_feature_filtering(uint16_t Nprocess, std::vector<std::string> files_vec, KmerGO_parameters p)
{
	std::string line, kmer;
	std::string input_file_name = files_vec[Nprocess];
	std::string onput_file_name_l = joinPath(p.features_path, "categorical_l_" + std::to_string(Nprocess) + ".txt");
	std::string onput_file_name_n = joinPath(p.features_path, "categorical_n_" + std::to_string(Nprocess) + ".txt");
	std::ifstream input_file(input_file_name);
	std::ofstream output_file_l(onput_file_name_l);
	std::ofstream output_file_n(onput_file_name_n);
	std::vector<short> sample_label;	// -1 inexistence, 0 groupA, 1 groupB
	std::string::size_type position;

	/* read head text and set labels */
	std::getline(input_file, line);
	mtx.lock();
	progress += line.size() + 1;	// update progress
	mtx.unlock();
	output_file_l << line << "\tASS\tLabel\n";	// write logical file head text
	output_file_n << line << "\tASS-l\tP\tASS-n\tLabel\n";	// write numerical file head text
	position = line.find("\t");
	line = line.substr(position + 1);	// remove string "k-mer"
	position = line.find("\t");
	while (position != line.npos)
	{
		std::string sample_name = line.substr(0, position);
		if (p.trait_information_map.find(sample_name) != p.trait_information_map.end())
		{
			if (p.trait_information_map[sample_name].compare(p.groupA_name) == 0) sample_label.push_back(0);
			else sample_label.push_back(1);
		}
		else sample_label.push_back(-1);
		line = line.substr(position + 1);
		position = line.find("\t");
	}
	if (p.trait_information_map.find(line) != p.trait_information_map.end())
	{
		if (p.trait_information_map[line].compare(p.groupA_name) == 0) sample_label.push_back(0);
		else sample_label.push_back(1);
	}
	else sample_label.push_back(-1);

	/* filtering */
	while (std::getline(input_file, line))
	{
		std::string src_line(line);
		position = line.find("\t");
		line = line.substr(position + 1);	// remove first k-mer
		position = line.find("\t");
		std::vector<float> group1, group2;
		int cnt = 0, tp = 0, tn = 0;
		char judge_label;
		while (position != line.npos)
		{
			std::string value_str = line.substr(0, position);
			float value_float = std::stof(value_str);
			/* add values to groups */
			if (sample_label[cnt] == 0)	// is groupA
			{
				group1.push_back(value_float);
				if (value_float > 0) ++tp;
			}
			else if(sample_label[cnt] == 1)	// is groupB
			{
				group2.push_back(value_float);
				if (value_float == 0) ++tn;
			}
			line = line.substr(position + 1);
			position = line.find("\t");
			++cnt;
		}
		/* calculation */
		float assL = ((float)tp / group1.size() + (float)tn / group2.size()) / 2;	// ass_l
		if (assL > 0.5) judge_label = 'A';
		else
		{
			judge_label = 'B';
			assL = 1 - assL;
		}
		if (assL >= p.assl)	// logical ASS
		{
			output_file_l << src_line << '\t' << std::to_string(assL) << '\t' << judge_label << '\n';
		}
		else	// numerical ASS
		{
			float pValue = ranksums(group1, group2);
			if (pValue < p.p_value)
			{
				std::vector<Logistic_Data> dataSet;
				for (float i : group1)
				{
					dataSet.push_back(Logistic_Data(std::vector<float>{i}, 0));	// 1 feature for group1
				}
				for (float i : group2)
				{
					dataSet.push_back(Logistic_Data(std::vector<float>{i}, 1));	// 1 feature for group2
				}
				std::vector<Logistic_Data> testdataSet(dataSet);
				Logistic model(dataSet);
				model.logisticRegression();
				model.predictClass(testdataSet);
				/* calculate the confusion matrix */
				tp = tn = 0;
				for (unsigned int i = 0; i < group1.size(); i++)	// caculate tp
				{
					if (testdataSet[i].cls == dataSet[i].cls) ++tp;
				}
				for (unsigned int i = 0; i < group1.size(); i++)	// caculate tn
				{
					if (testdataSet[i].cls == dataSet[i].cls) ++tn;
				}
				float assN = ((float)tp / group1.size() + (float)tn / group2.size()) / 2;	// ass_n
				if (assN >= p.assn)
				{
					output_file_n << src_line << '\t' << std::to_string(assL) << '\t' << pValue <<
						'\t' << assN << '\t' << judge_label << '\n';
				}
			}
		}
		mtx.lock();
		progress += src_line.size() + 1;	// update progress
		mtx.unlock();
	}
}

//--------------------------------------------------------------------------------------------------------------------
void continuous_feature_filtering(uint16_t Nprocess, std::vector<std::string> files_vec, KmerGO_parameters p)
{
	std::string line, kmer;
	std::string input_file_name = files_vec[Nprocess];
	std::string onput_file_name_l = joinPath(p.features_path, "continuous_l_" + std::to_string(Nprocess) + ".txt");
	std::string onput_file_name_n = joinPath(p.features_path, "continuous_n_" + std::to_string(Nprocess) + ".txt");
	std::ifstream input_file(input_file_name);
	std::ofstream output_file_l(onput_file_name_l);
	std::ofstream output_file_n(onput_file_name_n);
	std::vector<std::string> head_value;
	std::vector<bool> sample_valid;
	std::string::size_type position;

	/* read head text and set labels */
	std::getline(input_file, line);
	mtx.lock();
	progress += line.size() + 1;	// update progress
	mtx.unlock();
	output_file_l << line << "\tP\n";	// write logical file head text
	output_file_n << line << "\tP\tCorr\n";	// write numerical file head text
	position = line.find("\t");
	line = line.substr(position + 1);	// remove string "k-mer"
	position = line.find("\t");
	while (position != line.npos)
	{
		std::string sample_name = line.substr(0, position);
		head_value.push_back(sample_name);
		if (p.trait_information_map.find(sample_name) != p.trait_information_map.end()) sample_valid.push_back(true);
		else sample_valid.push_back(false);
		line = line.substr(position + 1);
		position = line.find("\t");
	}
	head_value.push_back(line);
	if (p.trait_information_map.find(line) != p.trait_information_map.end()) sample_valid.push_back(true);
	else sample_valid.push_back(false);

	/* filtering */
	while (std::getline(input_file, line))
	{
		std::string src_line(line);
		position = line.find("\t");
		line = line.substr(position + 1);	// remove first k-mer
		position = line.find("\t");
		std::vector<float> group1, group2, kmer_fre;
		int cnt = 0;
		while (position != line.npos)
		{
			std::string value_str = line.substr(0, position);
			float value_float = std::stof(value_str);
			/* logical filtering */
			if (sample_valid[cnt] == true)		// existence
			{
				kmer_fre.push_back(value_float);
				if (value_float == 0) group1.push_back(std::stof(p.trait_information_map[head_value[cnt]]));
				else group2.push_back(std::stof(p.trait_information_map[head_value[cnt]]));
			}
			line = line.substr(position + 1);
			position = line.find("\t");
			++cnt;
		}
		/* calculation */
		float pValue = ranksums(group1, group2);
		if (pValue < p.p_value)
		{
			output_file_l << src_line << '\t' << std::to_string(pValue) << '\n';
		}
		else
		{
			std::vector<float> group;
			group.insert(group.end(), group1.begin(), group1.end());
			group.insert(group.end(), group2.begin(), group2.end());
			float corr = spearman_correlation_coefficient(kmer_fre, group);
			if (corr >= p.corr_vaule)
			{
				output_file_n << src_line << '\t' << std::to_string(pValue) << '\t' << std::to_string(corr) << '\n';
			}
		}
		mtx.lock();
		progress += src_line.size() + 1;	// update progress
		mtx.unlock();
	}
}

//--------------------------------------------------------------------------------------------------------------------
short kmer_filtering(KmerGO_parameters& p, std::string& loginfo)
{
	uint64 Allfilesize = 0;
	std::vector<std::string> files_vec;

	getFiles(p.matrixes_path, files_vec, true);
	for (auto f : files_vec)
	{
		Allfilesize += getFileSize(f);
	}

	//start threads
	std::vector<std::thread> sub_thread;
	for (unsigned int i = 0; i < files_vec.size(); i++)
	{
		if (p.mode == 0)
		{
			sub_thread.push_back(std::thread(categorical_feature_filtering, i, files_vec, p));
		}
		else
		{
			sub_thread.push_back(std::thread(continuous_feature_filtering, i, files_vec, p));
		}
	}
	sub_thread.push_back(std::thread(show_Progress, Allfilesize));
	for (unsigned int i = 0; i < files_vec.size(); i++)
	{
		sub_thread[i].join();
	}
	sub_thread[files_vec.size()].join();

	switch (filtering_err_number)
	{
	default:
		break;
	}

	return 0;
}
