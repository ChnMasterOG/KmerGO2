/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2020-07-27
*/

#include "KmerGO.h"
#include "KmerGO_io.h"

static volatile uint64 progress = 0;
static volatile short union_err_number = 0;
static std::mutex mtx;

//--------------------------------------------------------------------------------------------------------------------
static void show_Progress(uint64 allkmernumber)
{
	printf("STEP k-mer union start!\n");
	while (progress < allkmernumber && union_err_number == 0)
	{
		printf("\rTotal: %.2f%%  ", (double)progress * 100 / allkmernumber);
	}
	if (union_err_number == 0)
	{
		printf("\rTotal: 100.00%%\n");
	}
}

//--------------------------------------------------------------------------------------------------------------------
void get_Son_Matrix(uint16_t Nprocess, KmerGO_parameters p, std::vector<std::string> kmerdb_list, 
	std::vector<std::string> beacon_kmer, std::vector<float> ALLfcounter)
{
	uint32 i;
	uint32 files_number = kmerdb_list.size();
	uint32 _kmer_length;

	CKMCFile* kmer_data_base = new CKMCFile[files_number];
	FILE* out_file;
	std::string out_file_name = joinPath(p.matrixes_path, "son_matrix_" + std::to_string(Nprocess) + ".txt");

	/* open k-mer data base */
	for (i = 0; i < files_number; i++)
	{
		if (!kmer_data_base[i].OpenForListing(joinPath(p.kmercountings_path, kmerdb_list[i])))
		{
			mtx.lock();
			union_err_number = -11;
			mtx.unlock();
			std::cout << "Thread-" << std::to_string(Nprocess) << "[Error]: Can not open some k-mer database in at least one thread." << std::endl;
			return;
		}
	}

	/* open the output file */
	if ((out_file = fopen(out_file_name.c_str(), "wb")) == NULL)
	{
		mtx.lock();
		union_err_number = -12;
		mtx.unlock();
		std::cout << "Thread-" << std::to_string(Nprocess) << "[Error]: Can not open output files in at least one thread." << std::endl;
		return;
	}
	setvbuf(out_file, NULL, _IOFBF, 1 << 24);

	/* get information from the k-mer database */
	_kmer_length = kmer_data_base[0].KmerLength();
	CKmerAPI kmer_object(_kmer_length);
	CKmerAPI beacon_start_kmer_object(_kmer_length);
	CKmerAPI beacon_end_kmer_object(_kmer_length);
	CKmerAPI minimum_kmer_object(_kmer_length);
	CKmerAPI maximum_kmer_object(_kmer_length);
	uint64 counter;

	/* create beacon k-mer object: start-kmer index is 0; end-kmer index is 1 */
	beacon_start_kmer_object.from_string(beacon_kmer[0]);
	beacon_end_kmer_object.from_string(beacon_kmer[1]);

	/* create minimum/maximum k-mer object */
	std::string temp_str = "";
	for (i = 0; i < _kmer_length - 1; i++) temp_str += "A";
	minimum_kmer_object.from_string(temp_str);	// is less than "A"*_kmer_length (modified kmer_api "operator<" function)
	temp_str = "";
	for (i = 0; i < _kmer_length; i++) temp_str += "T";
	maximum_kmer_object.from_string(temp_str);

	/* write head text */
	fwrite("k-mer", 1, 5, out_file);
	for (i = 0; i < files_number; i++)
	{
		fwrite("\t", 1, 1, out_file);
		fwrite(kmerdb_list[i].c_str(), 1, kmerdb_list[i].length(), out_file);
	}
	fwrite("\n", 1, 1, out_file);
	
	/* jump to corresponding k-mer and build k-mer loser tree */
	KmerLoserTree kmer_losertree(files_number, _kmer_length);
	uint32 end_number = 0;
	bool start_flag;
	for (i = 0; i < files_number; i++)
	{
		start_flag = false;
		while (kmer_data_base[i].ReadNextKmer(kmer_object, counter))
		{
			if (!kmer_object.operator<(beacon_start_kmer_object))
			{
				kmer_losertree.kmer_list[i].operator=(kmer_object);
				kmer_losertree.counter_list[i] = counter;
				start_flag = true;
				break;
			}
		}
		if (!start_flag || !kmer_object.operator<(beacon_end_kmer_object))	// no next k-mer
		{
			end_number++;
			kmer_losertree.kmer_list[i].operator=(maximum_kmer_object);
			kmer_losertree.counter_list[i] = 0;
		}
		else
		{
			mtx.lock();
			progress++;
			mtx.unlock();
		}
	}
	if (end_number == files_number)	// end of running
	{
		for (i = 0; i < files_number; i++) kmer_data_base[i].Close();
		fclose(out_file);
		delete[] kmer_data_base;
		return;
	}
	kmer_losertree.kmer_list[files_number].operator=(minimum_kmer_object);
	kmer_losertree.build();

	/* prepare variables for main loop */
	uint32 min_index;
	uint32 no_zero_counter[2] = { 0, 0 };
	uint32 no_zero_counter_thr[2] = { 1, 1 };
	std::vector<std::string> zero_matrix;
	std::vector<std::string> write_matrix;
	CKmerAPI min_kmer(_kmer_length);
	char temp_cstr[1024];

	if (p.single_step == false)
	{
		no_zero_counter_thr[0] = p.groupA_number * 0.2;
		if (!p.groupA_name.compare("")) no_zero_counter_thr[1] = 1;
		else no_zero_counter_thr[1] = (files_number - p.groupA_number) * 0.2;
	}

	for (i = 0; i < files_number; i++)
	{
		zero_matrix.push_back("\t0");
		write_matrix.push_back("\t0");
	}

	/* main loop */
	while (true)
	{
		min_index = kmer_losertree.getMIN();
		sprintf(temp_cstr, "\t%.4f", (double)kmer_losertree.counter_list[min_index] / ALLfcounter[min_index]);	// normalization
		if (min_kmer.operator==(kmer_losertree.kmer_list[min_index]))	// new k-mer is equal to last k-mer
		{
			write_matrix[min_index] = temp_cstr;
		}
		else
		{
			if (no_zero_counter[0] >= no_zero_counter_thr[0] || no_zero_counter[1] >= no_zero_counter_thr[1])	// write last line
			{
				min_kmer.to_string(temp_str);
				fwrite(temp_str.c_str(), 1, _kmer_length, out_file);
				for (i = 0; i < files_number; i++) fwrite(write_matrix[i].c_str(), 1, write_matrix[i].length(), out_file);
				fwrite("\n", 1, 1, out_file);
			}
			no_zero_counter[0] = no_zero_counter[1] = 0;
			write_matrix = zero_matrix;
			write_matrix[min_index] = temp_cstr;
		}
		if (p.single_step == false)	// count no zero numbers
		{
			if (!p.groupA_name.compare("")) no_zero_counter[0]++;
			else if(!p.trait_information_map[kmerdb_list[min_index]].compare(p.groupA_name)) no_zero_counter[0]++;
			else no_zero_counter[1]++;
		}
		else no_zero_counter[0]++;
		min_kmer.operator=(kmer_losertree.kmer_list[min_index]);
		if (!kmer_data_base[min_index].ReadNextKmer(kmer_object, counter) || !kmer_object.operator<(beacon_end_kmer_object))	// no next k-mer
		{
			end_number++;
			kmer_object = maximum_kmer_object;
			counter = 0;
		}
		else
		{
			mtx.lock();
			progress++;	// update progress
			mtx.unlock();
		}
		if (end_number == files_number)	// end of running
		{
			if (no_zero_counter[0] >= no_zero_counter_thr[0] || no_zero_counter[1] >= no_zero_counter_thr[1])	// write last line
			{
				min_kmer.to_string(temp_str);
				fwrite(temp_str.c_str(), 1, _kmer_length, out_file);
				for (i = 0; i < files_number; i++) fwrite(write_matrix[i].c_str(), 1, write_matrix[i].length(), out_file);
				fwrite("\n", 1, 1, out_file);
			}
			break;
		}
		kmer_losertree.kmer_list[min_index].operator=(kmer_object);
		kmer_losertree.counter_list[min_index] = counter;
		kmer_losertree.adjust(min_index);
	}

	/* close files */
	for (i = 0; i < files_number; i++) kmer_data_base[i].Close();
	fclose(out_file);

	delete[] kmer_data_base;

	return;
}

//--------------------------------------------------------------------------------------------------------------------
char num2ACGT(unsigned char num)
{
	if (num == 0) return 'A';
	else if (num == 1) return 'C';
	else if (num == 2) return 'G';
	else return 'T';
}

//--------------------------------------------------------------------------------------------------------------------
void distribute_kmer_nonuniformly(uint8_t step, uint32_t num, uint32_t maxnum, std::string& out_kmer)
{
	if (step == 4) return;
	if (num < 4 * maxnum / 10) {
		out_kmer += "A"; distribute_kmer_nonuniformly(step + 1, num, 4 * maxnum / 10, out_kmer);
	}
	else if (num < 7 * maxnum / 10) {
		out_kmer += "C"; distribute_kmer_nonuniformly(step + 1, num - 4 * maxnum / 10, 3 * maxnum / 10, out_kmer);
	}
	else if (num < 9 * maxnum / 10) {
		out_kmer += "G"; distribute_kmer_nonuniformly(step + 1, num - 7 * maxnum / 10, 2 * maxnum / 10, out_kmer);
	}
	else {
		out_kmer += "T"; distribute_kmer_nonuniformly(step + 1, num - 9 * maxnum / 10, maxnum / 10, out_kmer);
	}
}

//--------------------------------------------------------------------------------------------------------------------
void distribute_kmer_uniformly(uint8_t step, uint32_t num, std::string& out_kmer)
{
	if (step == 4) return;
	if (num % 4 == 0) out_kmer = "A" + out_kmer;
	else if (num % 4 == 1) out_kmer = "C" + out_kmer;
	else if (num % 4 == 2) out_kmer = "G" + out_kmer;
	else out_kmer = "T" + out_kmer;
	distribute_kmer_uniformly(step + 1, num / 4, out_kmer);
}

//--------------------------------------------------------------------------------------------------------------------
short kmer_union(KmerGO_parameters& p, std::string &loginfo)
{
	std::vector<std::string> files, ab_files, kmerdb_list, beacon_list, beacon_kmer;
	std::vector<uint64> ALLcounter;
	std::vector<float> ALLfcounter;
	uint32 files_number;
	uint32 i, j;
	uint64 k, Allkmernumber = 0;

	getFiles(p.kmercountings_path, ab_files, true);
	getFiles(p.kmercountings_path, files, false);
	files_number = files.size();

	for (i = 0; i < files_number; i++)
	{
		if (files[i].compare(files[i].length() - 8, 8, ".kmc_suf", 0, 8) == 0)
		{
			kmerdb_list.push_back(files[i].substr(0, files[i].length() - 8));
		}
	}

	if (kmerdb_list.size() * 2 < files_number)
	{
		loginfo = "The number of k-mer database is error.";
		return -1;
	}

	//count all k-mer counters and k-mer number
	CKMCFile kmer_data_base;
	CKmerAPI* kmer_object;
	uint64 temp, counter;
	uint32 _kmer_length = 4;
	for (i = 0; i < kmerdb_list.size(); i++)
	{
		temp = 0;
		if (!kmer_data_base.OpenForListing(joinPath(p.kmercountings_path, kmerdb_list[i])))
		{
			loginfo = "Can not open some k-mer database.";
			return -2;
		}
		_kmer_length = kmer_data_base.KmerLength();
		kmer_object = new CKmerAPI(_kmer_length);
		while (kmer_data_base.ReadNextKmer(*kmer_object, counter))
		{
			temp += counter;
			Allkmernumber++;
		}
		delete kmer_object;
		ALLcounter.push_back(temp);
		kmer_data_base.Close();
	}

	//divided to process_number blocks
	std::string temp_kmer;
	for (i = 0; i < p.process_number; i++)
	{
		temp_kmer = "";
		if (p.process_number < 40)	//distribute k-mer nonuniformly
			distribute_kmer_nonuniformly(0, i * 10000 / p.process_number, 10000, temp_kmer);
		else	//distribute k-mer uniformly
			distribute_kmer_uniformly(0, i * 256 / p.process_number, temp_kmer);
		for (j = 0; j < _kmer_length - 4; j++) temp_kmer += "A";
		beacon_list.push_back(temp_kmer);
	}
	temp_kmer = "TTTT";
	for (i = 0; i < _kmer_length - 4; i++) temp_kmer += "T";
	beacon_list.push_back(temp_kmer);
	
	// normalization coefficient
	short power = -5;
	for (i = 0; i < kmerdb_list.size(); i++)
	{
		j = 0;
		k = ALLcounter[i];
		while (k != 0)
		{
			k /= 10;
			j++;
		}
		if (power < (short)j - 4) power = j - 4;
	}
	for (i = 0; i < kmerdb_list.size(); i++)
	{
		if (power > 0) ALLfcounter.push_back(ALLcounter[i] / pow(10, power));
		else ALLfcounter.push_back(ALLcounter[i] * pow(10, -power));
	}

	//start threads
	std::vector<std::thread> sub_thread;
	for (i = 0; i < p.process_number; i++)
	{
		beacon_kmer.push_back(beacon_list[i]);
		beacon_kmer.push_back(beacon_list[i + 1]);
		sub_thread.push_back(std::thread(get_Son_Matrix, i, p, kmerdb_list, beacon_kmer, ALLfcounter));
		beacon_kmer.clear();
	}
	sub_thread.push_back(std::thread(show_Progress, Allkmernumber));
	for (i = 0; i < p.process_number; i++)
	{
		sub_thread[i].join();
	}
	sub_thread[p.process_number].join();

	switch (union_err_number)
	{
	case -11:
		loginfo = "Can not open some k-mer database in at least one thread.";
	case -12:
		loginfo = "Can not open output files in at least one thread.";
	default:
		break;
	}

	return union_err_number;
}
