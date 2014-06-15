#ifndef DEF_SAMPLE
#define DEF_SAMPLE

#include "config.h"
#include "def_transcript.h"


class trans_fragment
{
public:
	long get_trans_index();
	long get_copy_index();
	unsigned long get_pos_start();
	unsigned long get_pos_end();
	unsigned long get_length();
	trans_fragment* get_next();
	void assign_next(trans_fragment *next_frag);
	bool is_selected();

	trans_fragment* split_frag(unsigned long cut_pos); //split this fragment into two at the cutting position, return the new fragment
	bool size_selected(int window_low, int window_high, double prob_selection);

	trans_fragment(long trans_index, long copy_index, unsigned long pos_start, unsigned long pos_end);
	~trans_fragment();
private:
	long trans_index;
	long copy_index;
	unsigned long pos_start; //position in the original transcript
	unsigned long pos_end; //position in the original transcript
	unsigned long length;
	bool selected;
	trans_fragment *next;
};

class pool_fragment
{
public:
	void add_fragment(trans_fragment *new_frag);

	void generate_orig_pool(transcriptome *orig_transcriptome); //given a transcriptome (transcripts and number of copies), generate the original pool
	double random_fragmentation(long num_cut); //randomly cut the transcript fragments in the pool into smaller fragments
	long size_selection(int window_low, int window_high); //select the fragments with size in the window [low, high]

	void output_stat(string output_folder, int rnd_prefix, transcriptome *orig_transcriptome);

	pool_fragment();
	~pool_fragment();
private:
	vector <trans_fragment*> fragment_array;
};







#endif

