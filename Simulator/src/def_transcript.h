#ifndef DEF_TRANSCRIPT
#define DEF_TRANSCRIPT

#include "config.h"



class transcript
{
public:
	string get_name();
	long trans_index;
	unsigned long get_length();
	unsigned long get_num_copy();
	string get_chr();
	transcript* get_next();

	unsigned long random_num_copy();

	transcript(string name, string chr, unsigned long pos_start, unsigned long pos_end, unsigned long length, string strand, string geneNm);
	~transcript();
private:
	string name;
	unsigned long length;
	unsigned long num_copy;
	string chromosome;
	unsigned long pos_start;
	unsigned long pos_end;
	string strand;
	string geneNm;
	transcript *next;
};


class transcriptome
{
public:
	vector <transcript*> transcript_list;

	void input_transcripts_gaf(string gaf_filename, string output_folder);

	transcriptome();
	~transcriptome();
private:
	unsigned long total_num_copy;
};



#endif

