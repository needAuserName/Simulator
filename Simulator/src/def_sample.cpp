/*    
 *    def_sample.cpp		
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

#include "def_sample.h"
#include "general_functions.h"

string itostr(int t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}

trans_fragment::trans_fragment(long trans_index, long copy_index, unsigned long pos_start, unsigned long pos_end)
{
	this->trans_index = trans_index;
	this->copy_index = copy_index;
	this->pos_start = pos_start;
	this->pos_end = pos_end;
	this->length = pos_end - pos_start + 1;
	this->selected = false;
	this->next = NULL;
}

trans_fragment::~trans_fragment()
{

}

long trans_fragment::get_trans_index()
{
	return this->trans_index;
}

long trans_fragment::get_copy_index()
{
	return this->copy_index;
}

unsigned long trans_fragment::get_pos_start()
{
	return this->pos_start;
}

unsigned long trans_fragment::get_pos_end()
{
	return this->pos_end;
}

unsigned long trans_fragment::get_length()
{
	return this->length;
}

trans_fragment* trans_fragment::get_next()
{
	return this->next;
}

void trans_fragment::assign_next(trans_fragment *next_frag)
{
	this->next = next_frag;
}

bool trans_fragment::is_selected()
{
	return this->selected;
}

trans_fragment* trans_fragment::split_frag(unsigned long cut_pos)
{
	if (cut_pos <= 0)
		cut_pos = 1;
	if (cut_pos > length)
		cut_pos = length;

	unsigned long coordinate_cut = pos_start + cut_pos;

	//create a new fragment
	trans_fragment *new_fragment;
	if (coordinate_cut <= pos_end) 
		new_fragment = new trans_fragment (trans_index, copy_index, coordinate_cut, pos_end);
	else 
		new_fragment = NULL;
	
	//modify this fragment
	pos_end = coordinate_cut-1;
	length = pos_end - pos_start + 1;
//	this->next = new_fragment;

	return new_fragment;
}

bool trans_fragment::size_selected(int window_low, int window_high, double prob_selection)
{
	if (this->length >= window_low && this->length <= window_high)
	{
		uniform_real_distribution<double> distribution_unif(0.0,1.0);

		double rand_num_unif = distribution_unif(rand_generator);

		if (rand_num_unif <= prob_selection)
		{
			// selected
			selected = true;
			return true;
		} 
		else
		{
			selected = false;
			return false;
		}
		
	} 
	else
	{
		selected = false;
		return false;
	}
}

pool_fragment::pool_fragment()
{
	fragment_array.reserve(const_num_cut_per_round * const_max_sequencing_round);
}

pool_fragment::~pool_fragment()
{
	for (unsigned long tmp_loop = 0; tmp_loop < fragment_array.size(); ++tmp_loop)
	{
		delete fragment_array[tmp_loop];
	}
	fragment_array.clear();
}

void pool_fragment::add_fragment(trans_fragment *new_frag)
{
	if (!new_frag)
		return;

	fragment_array.push_back(new_frag);

	return;
}

void pool_fragment::generate_orig_pool(transcriptome *orig_transcriptome)
{
	transcript *cur_trans;
	trans_fragment *new_frag;
	
	for (unsigned long trans_loop = 0; trans_loop < orig_transcriptome->transcript_list.size(); ++trans_loop)
	{
		cur_trans = orig_transcriptome->transcript_list[trans_loop];

		for (int frag_loop = 1; frag_loop <= cur_trans->get_num_copy(); ++frag_loop)
		{
			new_frag = new trans_fragment(trans_loop, frag_loop, 1, cur_trans->get_length());
			add_fragment(new_frag);
		}
	}

	return;
}

double pool_fragment::random_fragmentation(long num_cut)
{
	if (fragment_array.empty())
	{
		return 0;
	}

	uniform_real_distribution<double> distribution(0.0,1.0);
	double rand_num_unif;
	unsigned long last_pos = 0, rand_offset, rand_cut_position;
	trans_fragment *cut_fragment;

	// new: probability of cutting 
	// fit a 2nd degree polynomial for the cdf, that passes (0, 0) and (target_low*sigma + target_high, 1), y = c x^2
	// probability = c x^2, x is the length
	// cut a fragment with probability of 1 if it's longer than target_low*sigma + target_high
	const double const_sigma = 6;
	double coefficient_c = 1.0 / pow(double(const_fragmentation_target_high + const_fragmentation_target_low * const_sigma), 2);


	//randomly cut
	for (long cut_loop = 0, try_cnt = 0; cut_loop < num_cut && try_cnt < num_cut * 5; ++try_cnt)
	{
		//step 1, randomly choose the fragment to be cut
		rand_num_unif = distribution(rand_generator);
		rand_offset = floor(fragment_array.size() * rand_num_unif);

		//last_pos = (last_pos + rand_offset) % fragment_array.size();
		last_pos = rand_offset % fragment_array.size();
		cut_fragment = fragment_array[last_pos];

		//new on 8/8/2013: determine to cut or not based on the length of the fragment
		bool flag_cutting = false;
//		if (cut_fragment->get_length() <= const_fragmentation_target_low)
//		{
//			flag_cutting = false; // do not cut a short fragment
//		}
		if (cut_fragment->get_length() >= const_fragmentation_target_high + const_fragmentation_target_low * const_sigma)
		{
			flag_cutting = true; // cut a long fragment
		}
		else
		{			
			rand_num_unif = distribution(rand_generator);
			if (rand_num_unif <= coefficient_c * pow(double(cut_fragment->get_length()), 2))
			{
				flag_cutting = true; // bad luck, will cut you
			} 
			else
			{
				flag_cutting = false;
			}
		}

		if (flag_cutting)
		{
			++cut_loop;

			//step 2, randomly pick the cutting position
			rand_num_unif = distribution(rand_generator);
			rand_cut_position = floor(0.5 + cut_fragment->get_length() * rand_num_unif);

			//step 3, modify this fragment and add the new fragment
			add_fragment(cut_fragment->split_frag(rand_cut_position));
		}
	}

	long cnt_short_fragment = 0, cnt_in_target_fragment = 0, cnt_long_fragment = 0;
	for (unsigned long loop_fragment_array = 0; loop_fragment_array < fragment_array.size(); ++ loop_fragment_array)
	{
		if (fragment_array[loop_fragment_array]->get_length() < const_fragmentation_target_low)
		{
			++cnt_short_fragment;
		}
		else if (fragment_array[loop_fragment_array]->get_length() >= const_fragmentation_target_low && fragment_array[loop_fragment_array]->get_length() <= const_fragmentation_target_high)
		{
			++cnt_in_target_fragment;
		}
		else if (fragment_array[loop_fragment_array]->get_length() > const_fragmentation_target_high)
		{
			++cnt_long_fragment;
		}
	}
	cout << "# of all fragments = " << fragment_array.size() << "\t# of short fragments = " << cnt_short_fragment << "\t# of fragments in targeted window = " << cnt_in_target_fragment
		<< "\t# of long fragments = " << cnt_long_fragment << "\tratio in window = " << double(cnt_in_target_fragment) / fragment_array.size() << endl;
	
	return double(cnt_short_fragment + cnt_in_target_fragment) / fragment_array.size(); //return the ratio of fragments that should not be fragmented any more
}


long pool_fragment::size_selection(int window_low, int window_high)
{
	if (fragment_array.empty())
	{
		return 0;
	}

	//count the size of the library
	long num_frag_in_window = 0;
	for (unsigned long frag_loop = 0; frag_loop < fragment_array.size(); ++frag_loop)
	{
		if (fragment_array[frag_loop]->get_length() >= window_low && fragment_array[frag_loop]->get_length() <= window_high)
		{
			++num_frag_in_window;
		}
	}
	cout << "library size (total # of fragments in the size-selection window) = " << num_frag_in_window << endl;
	double prob_selection = double(const_desired_num_reads) / num_frag_in_window;

	long num_selected_frag = 0;
	for (unsigned long frag_loop = 0; frag_loop < fragment_array.size(); ++frag_loop)
	{
		if (fragment_array[frag_loop]->size_selected(window_low, window_high, prob_selection))
		{
			++num_selected_frag;
		}
	}

	cout << "# of sampled fragments = " << num_selected_frag << endl;

	return num_selected_frag;
}

void pool_fragment::output_stat(string output_folder, int rnd_prefix, transcriptome *orig_transcriptome)
{
	string filename;
	ofstream stat_file, length_file_selected, length_file_not_selected, length_file_total, frag_file_selected, trans_file;
	unsigned long num_selected = 0, num_not_selected = 0, length;

	filename = output_folder + itostr(rnd_prefix) + "_length_total.txt"; length_file_total.open(filename.c_str());
	long *fraglength = new long[max_translength + 1];
	for (unsigned int loop = 0; loop <= max_translength; loop++)
	{
		fraglength[loop] = 0;
	}
	vector <trans_fragment*> sequencedFrag_array;
	for (unsigned long frag_loop = 0; frag_loop < fragment_array.size(); ++frag_loop)
	{
		trans_fragment *cur_fragment = fragment_array[frag_loop];
		if (cur_fragment->is_selected())
		{
			sequencedFrag_array.push_back(cur_fragment);
		}
		fraglength[cur_fragment->get_length()]++;
	}
	for (unsigned int loop = 1; loop <= max_translength; loop++)
	{
		length_file_total << loop << "\t" << fraglength[loop] << endl;
	}
	delete [] fraglength;

	// sort the selected sequenced fragments
	void** sortlist_frags = new void* [sequencedFrag_array.size() + 2];
	double* sortkey_frags = new double [sequencedFrag_array.size() + 2];
	for (unsigned long frag_loop = 0; frag_loop < sequencedFrag_array.size(); ++frag_loop)
	{
		sortlist_frags[frag_loop + 1] = (void*)sequencedFrag_array[frag_loop];
		sortkey_frags[frag_loop + 1] = sequencedFrag_array[frag_loop]->get_pos_start();
	}
	mergeSort_general(sortlist_frags, sortkey_frags, sequencedFrag_array.size());
	for (unsigned long frag_loop = 0; frag_loop < sequencedFrag_array.size(); ++frag_loop)
	{
		sequencedFrag_array[frag_loop] = (trans_fragment*) sortlist_frags[frag_loop + 1];
	}
	for (unsigned long frag_loop = 0; frag_loop < sequencedFrag_array.size(); ++frag_loop)
	{
		sortlist_frags[frag_loop + 1] = (void*)sequencedFrag_array[frag_loop];
		sortkey_frags[frag_loop + 1] = sequencedFrag_array[frag_loop]->get_copy_index();
	}
	mergeSort_general(sortlist_frags, sortkey_frags, sequencedFrag_array.size());
	for (unsigned long frag_loop = 0; frag_loop < sequencedFrag_array.size(); ++frag_loop)
	{
		sequencedFrag_array[frag_loop] = (trans_fragment*) sortlist_frags[frag_loop + 1];
	}
	for (unsigned long frag_loop = 0; frag_loop < sequencedFrag_array.size(); ++frag_loop)
	{
		sortlist_frags[frag_loop + 1] = (void*)sequencedFrag_array[frag_loop];
		sortkey_frags[frag_loop + 1] = sequencedFrag_array[frag_loop]->get_trans_index();
	}
	mergeSort_general(sortlist_frags, sortkey_frags, sequencedFrag_array.size());
	for (unsigned long frag_loop = 0; frag_loop < sequencedFrag_array.size(); ++frag_loop)
	{
		sequencedFrag_array[frag_loop] = (trans_fragment*) sortlist_frags[frag_loop + 1];
	}
	delete [] sortlist_frags;
	delete [] sortkey_frags;

	// output fragment length stats
	filename = output_folder + "stat.txt"; stat_file.open(filename.c_str(), fstream::app);
	filename = output_folder + itostr(rnd_prefix) + "_length_selected.txt"; length_file_selected.open(filename.c_str());
	filename = output_folder + itostr(rnd_prefix) + "_length_not_selected.txt"; length_file_not_selected.open(filename.c_str());
	filename = output_folder + itostr(rnd_prefix) + "_fraginfo_selected.txt"; frag_file_selected.open(filename.c_str());

	unsigned long frag_loop = 0, frag_prev_loop, frag_post_loop;
	while(frag_loop < sequencedFrag_array.size())
	{
		if (sequencedFrag_array[frag_loop]->get_pos_start() > 1)
		{
			length = sequencedFrag_array[frag_loop]->get_pos_start() - 1;
			length_file_not_selected << length << endl;
			++num_not_selected;
		}
		length_file_selected << sequencedFrag_array[frag_loop]->get_length() << endl;
		frag_file_selected << sequencedFrag_array[frag_loop]->get_trans_index() << "\t" << sequencedFrag_array[frag_loop]->get_copy_index() << "\t" << sequencedFrag_array[frag_loop]->get_pos_start() << "\t" << sequencedFrag_array[frag_loop]->get_pos_end() << "\t" << sequencedFrag_array[frag_loop]->get_length() << endl;
		++num_selected;

		if (frag_loop == sequencedFrag_array.size() - 1)
		{
			if (sequencedFrag_array[frag_loop]->get_pos_end() < orig_transcriptome->transcript_list[sequencedFrag_array[frag_loop]->get_trans_index()]->get_length())
			{
				length = orig_transcriptome->transcript_list[sequencedFrag_array[frag_loop]->get_trans_index()]->get_length() - sequencedFrag_array[frag_loop]->get_pos_end();
				length_file_not_selected << length << endl;
				++num_not_selected;
			}
			break;
		}

		frag_prev_loop = frag_loop;
		frag_post_loop = frag_prev_loop + 1;
		while(frag_post_loop < sequencedFrag_array.size())
		{
			if (sequencedFrag_array[frag_post_loop]->get_trans_index() != sequencedFrag_array[frag_prev_loop]->get_trans_index())
			{
				if (sequencedFrag_array[frag_prev_loop]->get_pos_end() < orig_transcriptome->transcript_list[sequencedFrag_array[frag_prev_loop]->get_trans_index()]->get_length())
				{
					length = orig_transcriptome->transcript_list[sequencedFrag_array[frag_prev_loop]->get_trans_index()]->get_length() - sequencedFrag_array[frag_prev_loop]->get_pos_end();
					length_file_not_selected << length << endl;
					++num_not_selected;
				}
				frag_loop = frag_post_loop;
				break;
			}
			else
			{
				if (sequencedFrag_array[frag_post_loop]->get_copy_index() != sequencedFrag_array[frag_prev_loop]->get_copy_index())
				{
					if (sequencedFrag_array[frag_prev_loop]->get_pos_end() < orig_transcriptome->transcript_list[sequencedFrag_array[frag_prev_loop]->get_trans_index()]->get_length())
					{
						length = orig_transcriptome->transcript_list[sequencedFrag_array[frag_prev_loop]->get_trans_index()]->get_length() - sequencedFrag_array[frag_prev_loop]->get_pos_end();
						length_file_not_selected << length << endl;
						++num_not_selected;
					}
					frag_loop = frag_post_loop;
					break;
				}
				else
				{
					if (sequencedFrag_array[frag_post_loop]->get_pos_start() > (sequencedFrag_array[frag_prev_loop]->get_pos_end() + 1))
					{
						length = sequencedFrag_array[frag_post_loop]->get_pos_start() - sequencedFrag_array[frag_prev_loop]->get_pos_end() - 1;
						length_file_not_selected << length << endl;
						++num_not_selected;
					}
				}
			}

			length_file_selected << sequencedFrag_array[frag_post_loop]->get_length() << endl;
			frag_file_selected << sequencedFrag_array[frag_post_loop]->get_trans_index() << "\t" << sequencedFrag_array[frag_post_loop]->get_copy_index() << "\t" << sequencedFrag_array[frag_post_loop]->get_pos_start() << "\t" << sequencedFrag_array[frag_post_loop]->get_pos_end() << "\t" << sequencedFrag_array[frag_post_loop]->get_length() << endl;
			++num_selected;
			frag_loop = frag_post_loop + 1;
			frag_prev_loop = frag_post_loop;
			frag_post_loop++;
		}
		
	}


	stat_file << rnd_prefix << "\t" << num_selected + num_not_selected << "\t" << num_selected << "\t" << num_not_selected << endl;

// 	num_selected = num_not_selected = 0;
// 	for (frag_loop = 0; frag_loop < fragment_array.size(); ++frag_loop)
// 	{
// 		trans_fragment *cur_fragment = fragment_array[frag_loop];
// 		if (cur_fragment->is_selected())
// 		{
// 			++num_selected;
// 		}
// 		else
// 		{
// 			++num_not_selected;
// 		}
// 	}
// 	stat_file << rnd_prefix << "\t" << num_selected + num_not_selected << "\t" << num_selected << "\t" << num_not_selected << endl;

	// output stats showing whether each copy has been sampled
// 	filename = output_folder + itostr(rnd_prefix) + "_trans_info.txt"; trans_file.open(filename.c_str());
// 	bool sampled = false;
// 	for (unsigned long trans_loop = 0; trans_loop < orig_transcriptome->transcript_list.size(); trans_loop++)
// 	{
// 		for (int copy_loop = 1; copy_loop <= orig_transcriptome->transcript_list[trans_loop]->get_num_copy(); ++copy_loop)
// 		{
// 			sampled = false;
// 			trans_file << orig_transcriptome->transcript_list[trans_loop]->get_name() << "\t" << orig_transcriptome->transcript_list[trans_loop]->get_length() << "\t" << orig_transcriptome->transcript_list[trans_loop]->trans_index << "\t" << copy_loop << "\t";
// 			for (frag_loop = 0; frag_loop < sequencedFrag_array.size(); frag_loop++)
// 			{
// 				if (sequencedFrag_array[frag_loop]->get_trans_index() > trans_loop)
// 				{
// 					break;
// 				}
// 				if (sequencedFrag_array[frag_loop]->get_trans_index() == trans_loop && sequencedFrag_array[frag_loop]->get_copy_index() == copy_loop)
// 				{
// 					sampled = true;
// 					break;
// 				}
// 			}
// 			trans_file << sampled << endl;
// 		}
// 	}
// 	sequencedFrag_array.clear();
	
	stat_file.close();
	length_file_selected.close();
	length_file_not_selected.close();
	length_file_total.close();
	frag_file_selected.close();
//	trans_file.close();

	return;
}







