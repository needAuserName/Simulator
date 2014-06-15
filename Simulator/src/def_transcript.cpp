#include "def_transcript.h"

transcript::transcript(string name, string chr, unsigned long pos_start, unsigned long pos_end, unsigned long length, string strand, string geneNm)
{
	this->name = name;
	this->length = length;
	this->chromosome = chr;
	this->pos_start = pos_start;
	this->pos_end = pos_end;
	this->strand = strand;
	this->geneNm = geneNm;

	this->num_copy = 0;
	this->next = NULL;
}

transcript::~transcript()
{
	//nothing
}

string transcript::get_name()
{
	return this->name;
}


unsigned long transcript::get_length()
{
	return this->length;
}

unsigned long transcript::get_num_copy()
{
	return this->num_copy;
}

string transcript::get_chr()
{
	return this->chromosome;
}

transcript* transcript::get_next()
{
	return this->next;
}

//randomly assign the number of copies of this transcript in the transcriptome
unsigned long transcript::random_num_copy()
{
	//stage 1, determine whether this transcript is present

	uniform_real_distribution<double> distribution_unif(0.0,1.0);

	double rand_num_unif = distribution_unif(rand_generator);

	if (rand_num_unif <= const_prob_trans_present)
	{
		// present
	} 
	else
	{
		// not present
		this->num_copy = 0;
		return 0;
	}

	//stage 2, determine the number of copies for this transcript 

	negative_binomial_distribution<int> distribution_nb(const_num_copy_NB_r, const_num_copy_NB_p);

	int rand_num_nb = distribution_nb(rand_generator);
	if (rand_num_nb <= const_max_num_copy)
	{
		this->num_copy = rand_num_nb;
	}
	else
	{
		this->num_copy = floor( (const_num_copy_NB_r - 1.0) * (1.0 - const_num_copy_NB_p) / const_num_copy_NB_p ); //mode of this NB
	}

	return this->num_copy;
}



transcriptome::transcriptome()
{
	total_num_copy = 0;
}

transcriptome::~transcriptome()
{
	for (unsigned long tmp_loop = 0; tmp_loop < transcript_list.size(); ++tmp_loop)
	{
		delete transcript_list[tmp_loop];
	}
	transcript_list.clear();
}


void transcriptome::input_transcripts_gaf(string gaf_filename, string output_folder)
{
	string cur_exonboundary, tmpStr, trans_name, trans_chr, trans_strand, trans_geneNm;
	ifstream gaf_file;
	long this_boundary, last_boundary, trans_pos_start, trans_pos_end; 
	unsigned long trans_length;
	transcript *new_trans;
	bool validTrans = false;

	gaf_file.open(gaf_filename.c_str());
	max_translength = 0;

	while (getline(gaf_file, tmpStr, '\t'))
	{
		getline(gaf_file, trans_name, '\t');
		for (int tmpLoop = 1; tmpLoop <= 12; ++tmpLoop)
			getline(gaf_file, tmpStr, '\t');

		getline(gaf_file, trans_chr, ':');
		getline(gaf_file, cur_exonboundary, ':');
		getline(gaf_file, trans_strand, '\t');
		getline(gaf_file, trans_geneNm, '\t');
		getline(gaf_file, tmpStr);

		if (transcript_list.size() > 0 && trans_name.compare(transcript_list[transcript_list.size()-1]->get_name()) == 0)
		{
			continue;
		}

		trans_pos_start = MAX_NUMBER;
		trans_pos_end = 0;
		last_boundary = -1;
		trans_length = 0;
		while (cur_exonboundary.empty() == false)
		{
			size_t found = cur_exonboundary.find_first_of(" -,\n");
			if (found != string::npos)
			{
				tmpStr = cur_exonboundary.substr(0, found);
				cur_exonboundary.erase(0, found+1);
			}
			else
			{
				tmpStr = cur_exonboundary.substr(0);
				cur_exonboundary.clear();
			}

			this_boundary = atol(tmpStr.c_str());
			this_boundary = abs(this_boundary);

// 			if (new_annogene->exonBoundary.size() >= new_annogene->exonBoundary.capacity())
// 				new_annogene->exonBoundary.reserve(new_annogene->exonBoundary.capacity() + 10);
// 			new_annogene->exonBoundary.push_back(tmpboundary);

			if (last_boundary > 0)
			{
				assert(this_boundary >= last_boundary);
				trans_length += this_boundary - last_boundary + 1;
				last_boundary = -1;
			} 
			else
			{
				last_boundary = this_boundary;
			}

			if (this_boundary < trans_pos_start)
				trans_pos_start = this_boundary;
			if (this_boundary > trans_pos_end)
				trans_pos_end = this_boundary;
		}
		assert(last_boundary < 0);

		validTrans = false;
		if (trans_chr[0] == 'c' && trans_chr[1] == 'h' && trans_chr[2] == 'r' && (trans_chr[3] >= '0' && trans_chr[3] <= '9' || trans_chr[3] == 'X' || trans_chr[3] == 'Y'))
		{
			validTrans = true;
		}
		if (!trans_name.empty() && !trans_chr.empty() && trans_pos_start <= trans_pos_end && trans_length > 0 && validTrans == true)
		{
			new_trans = new transcript (trans_name, trans_chr, trans_pos_start, trans_pos_end, trans_length, trans_strand, trans_geneNm);
			if (transcript_list.size() >= transcript_list.capacity())
			{
				transcript_list.reserve(transcript_list.capacity() + 5000);
			}
			transcript_list.push_back(new_trans);
			new_trans->trans_index = transcript_list.size() - 1;
			total_num_copy += new_trans->random_num_copy();
		}
	}


	string filename = output_folder + "stat.txt";
	ofstream stat_file(filename.c_str(), fstream::app);
	stat_file << transcript_list.size() << "\t" << total_num_copy << endl;
	stat_file.close();

	filename = output_folder + "transcript.txt";
	ofstream transcript_file(filename.c_str());
	transcript_file << "transindex\tchromosome\tname\tlength\t#copy" << endl;
	for (unsigned long tmp_loop = 0; tmp_loop < transcript_list.size(); ++tmp_loop)
	{
		transcript_file << transcript_list[tmp_loop]->trans_index << "\t" << transcript_list[tmp_loop]->get_chr() << "\t" << transcript_list[tmp_loop]->get_name() << "\t" << transcript_list[tmp_loop]->get_length() << "\t" << transcript_list[tmp_loop]->get_num_copy() << endl;
		if (transcript_list[tmp_loop]->get_length() > max_translength)
		{
			max_translength = transcript_list[tmp_loop]->get_length();
		}
	}
	transcript_file.close();

	return;
}



