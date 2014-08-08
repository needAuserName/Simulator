/*    
 *    sim_gap.cpp		
 *
 *    Copyright (C) 2014 University of Kentucky and
 *                       Yan Huang
 *
 *    Authors: Yan Huang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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