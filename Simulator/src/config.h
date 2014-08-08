/*    
 *    config.cpp		
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

/////////////////////////////////////////////////////////////////////////////
// Simulate the fragmentation and size selection of an RNA-seq experiment

#ifndef CONFIG
#define CONFIG

//#define UNIX

#ifdef UNIX
#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <sstream>
#include <cmath>
#include <ctime>
#include <random>
#else
#include <fstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <assert.h>
#include <sstream>
#include <time.h>
#include <random>
#endif

using namespace std;


const long MAX_NUMBER = 2000000000;

//////////////////////////////////////////////////////////////////////////
// General parameters
const long const_desired_num_reads = 5e7; // expected output size, reflecting the sampling depth

//////////////////////////////////////////////////////////////////////////
// Parameters for random fragmentation
const int const_max_sequencing_round = 50;
const long const_num_cut_per_round = 1e7;
// const long const_num_cut_per_round = 5e2;
// const long const_desired_num_reads = 1e2;
const int const_fragmentation_target_low = 200;  // this window defines the expected length range of the mRNA fragments after fragmentation,
const int const_fragmentation_target_high = 500; // and is different from the size selection window below; use this window to measure the quality of fragmentation
const double const_thresh_fragment_in_target_ratio = 0.9; // end the fragmentation if this ratio of the total fragments have been in the targeted window


//////////////////////////////////////////////////////////////////////////
// Parameters for size selection
const int const_size_window_low = 150; 
const int const_size_window_high = 350;
//const double const_prob_seqfrag_selected = 0.7;


//////////////////////////////////////////////////////////////////////////
// Parameters for transcript expression
const double const_prob_trans_present = 0.9;
const double const_num_copy_NB_r = 10.0;
const double const_num_copy_NB_p = 0.05;
const int const_max_num_copy = 2000;
//mean: r (1-p) / p
//mode: floor( (r-1) (1-p) / p )
//variance: r (1-p) / p^2



extern default_random_engine rand_generator;

extern long max_translength;



#endif
