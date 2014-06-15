#include "config.h"
#include "def_transcript.h"
#include "def_sample.h"

default_random_engine rand_generator;


string output_folder = "result/";
long max_translength;


int main (int argc, char *argv[])
{	
// 	srand(time(NULL));
// 	rand_generator.seed (rand());
	rand_generator.seed (time(NULL));

	//prepare the transcriptome
	transcriptome custom_transtome;
	pool_fragment custom_pool;
	long num_fragment_selected, num_fragment_selected_prev;

	custom_transtome.input_transcripts_gaf("trans.gaf", output_folder);
	custom_pool.generate_orig_pool(&custom_transtome);

	custom_pool.output_stat(output_folder, 0, &custom_transtome);
	num_fragment_selected_prev = 0;
// 	for (int seq_round = 1; seq_round <= const_max_sequencing_round; ++seq_round)
// 	{
// 		custom_pool.random_fragmentation(const_num_cut_per_round);
// 		num_fragment_selected = custom_pool.size_selection(const_size_window_low, const_size_window_high);
// 		custom_pool.output_stat(output_folder, seq_round, &custom_transtome);
// 
// 		if (num_fragment_selected > const_desired_num_reads)
// 			break;
// 		if (num_fragment_selected < num_fragment_selected_prev)
// 			break;
// 		num_fragment_selected_prev = num_fragment_selected;
// 	}
	int seq_round;
	for (seq_round = 1; seq_round <= const_max_sequencing_round; ++seq_round)
	{
		cout << seq_round << "\t";
		double cur_ratio_frag_in_target = custom_pool.random_fragmentation(const_num_cut_per_round);
		if (cur_ratio_frag_in_target > const_thresh_fragment_in_target_ratio)
		{
			break;
		}
	}
	num_fragment_selected = custom_pool.size_selection(const_size_window_low, const_size_window_high);
	custom_pool.output_stat(output_folder, seq_round, &custom_transtome);


	string filename = output_folder + "stat.txt";
	ofstream stat_file(filename.c_str(), fstream::app);
	stat_file << "simulation finished" << endl << endl << endl;
	stat_file.close();

	return 0;
}