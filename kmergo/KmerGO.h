/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2020-07-27
*/

#ifndef _KEMRGO_H
#define _KMERGO_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include <math.h>
#include "../stats_ml/stats_ml.h"
#include "../kmc_api/kmc_file.h"

#define KMERGO_VER	"2.0.0"
#define KMERGO_DATE	"2022-05-09"

#define MAX_FILES	65535
#define MAX_THREAD	256
#define MAX_READ	1024

typedef enum
{
	NO_SELECTED = 0,
	KMC3,
	UNION,
	FILTERING,
	CAP3,
	ALL,
}KmerGO_tools;

typedef struct
{
	//input parameters
	unsigned short k_value = 40;
	unsigned int ci_value = 2;
	unsigned int cs_value = 65535;
	unsigned short process_number = 24;
	unsigned char mode = 0;
	float assl = 0.8;
	float p_value = 0.01;
	float assn = 0.8;
	float corr_vaule = 0.8;
	std::string samples_path = "./";
	std::string kmercountings_path = "./";
	std::string matrixes_path = "./";
	std::string features_path = "./";
	std::string results_path = "./";
	std::string temp_path = "./";
	std::string csv_path = "";
	//running information
	bool single_step = true;
	//group information
	std::string groupA_name;
	uint32 groupA_number;
	std::map<std::string, std::string> trait_information_map;
}KmerGO_parameters;

class KmerLoserTree
{
private:
	std::vector<uint32> ls;	//a winner is in ls[0], losers are in others
	void swap(uint32* a, uint32* b);
public:
	KmerLoserTree(uint32 c, uint32 k);
	~KmerLoserTree();
	uint32 capacity;
	std::vector<CKmerAPI> kmer_list;
	std::vector<uint32> counter_list;
	void adjust(uint32 s);
	void build();
	uint32 getMIN();
};

short kmer_counting(KmerGO_parameters& p, std::string& loginfo);
short kmer_union(KmerGO_parameters& p, std::string& loginfo);
short kmer_filtering(KmerGO_parameters& p, std::string& loginfo);
short kmer_assembly(KmerGO_parameters& p, std::string& loginfo);

bool Read_CSV(std::string csv_path, std::map<std::string, std::string>& stringMap);

#endif
