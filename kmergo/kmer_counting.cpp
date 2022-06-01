/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2020-07-27
*/

#include "KmerGO.h"
#include "KmerGO_io.h"

//--------------------------------------------------------------------------------------------------------------------
short kmer_counting(KmerGO_parameters& p, std::string& loginfo)
{
	std::vector<std::string> suffix{ ".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz", ".fastq", ".fq", ".fastq.gz", ".fq.gz" };
	std::vector<std::string> kmc_ftype{ "-fm", "-fm", "-fm", "-fm", "-fm", "-fm", "-fq", "-fq", "-fq", "-fq" };

	std::string current_path, cmd, file_type = "?";
	std::vector<std::string> kmer_files, ap_kmer_files;
	uint32 r, i, j;

	current_path = getCWD();
	getFiles(p.samples_path, kmer_files, false);
	getFiles(p.samples_path, ap_kmer_files, true);

	for (i = 0; i < kmer_files.size(); i++)
	{
		std::cout << kmer_files[i] << std::endl;
		for (j = 0; j < suffix.size(); j++)
		{
			if (!suffix[j].compare(0, suffix[j].length(), kmer_files[i], kmer_files[i].length() - suffix[j].length(), suffix[j].length()))
			{
				file_type = kmc_ftype[j];
				break;
			}
		}
		if (file_type.compare("?") == 0) continue;	// skip illegal files
		cmd = joinPath(current_path, "bin/kmc -k") + std::to_string(p.k_value) + " -ci" + std::to_string(p.ci_value) +
			" -cs" + std::to_string(p.cs_value) + " " + file_type + " " + ap_kmer_files[i] + " " + 
			joinPath(p.kmercountings_path, kmer_files[i]) + " " + p.temp_path;
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "Running kmc error.";
			return 2;
		}
		cmd = joinPath(current_path, "bin/kmc_tools transform ") + joinPath(p.kmercountings_path, kmer_files[i]) + " sort " +
			joinPath(p.kmercountings_path, kmer_files[i]) + "_sort";
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "Running kmc_tools error.";
			return 3;
		}
		cmd = "mv " + joinPath(p.kmercountings_path, kmer_files[i]) + "_sort.kmc_suf " + 
			joinPath(p.kmercountings_path, kmer_files[i]) + ".kmc_suf";
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "Moving k-mer database error.";
			return 4;
		}
		cmd = "mv " + joinPath(p.kmercountings_path, kmer_files[i]) + "_sort.kmc_pre " + 
			joinPath(p.kmercountings_path, kmer_files[i]) + ".kmc_pre";
		r = system(cmd.c_str());
		if (r)
		{
			loginfo = "Moving k-mer database error.";
			return 4;
		}
	}

	return 0;
}


