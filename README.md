Simulator
=========

simulate RNA-seq reads following RNA-seq protocal

===============

[Yan Huang](http://protocols.netlab.uky.edu/~yan) \(yan at netlab dot uky dot edu\)

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Compilation & Installation](#compilation)
* [Usage](#usage)
* [Example](#example)
* [Simulation](#simulation)
* [Generate Transcript-to-Gene-Map from Trinity Output](#gen_trinity)
* [Differential Expression Analysis](#de)
* [Acknowledgements](#acknowledgements)

* * *

## <a name="introduction"></a> Introduction

The simulator is designed to mimic a real RNA-seq experiment and 
generates fragments from provided transcript copies. The simulation
process consists of three steps: (1) Build a synthetic transcriptome 
by randomly assign copy numbers to all the genes and isoforms in 
the annotation database and set this as the true profiles. 
(2) Randomly cut the transcripts in the synthetic transcriptome into 
small fragments and dynamically check the lengths of the generated 
fragments. Fragments with lengths in a certain range (e.g. [150bp; 350bp]) 
are selected with probability to construct the sequencing library. 
This step stops when the number of fragments in the library exceeds
the pre-specied sequencing depth. (3) 2*75bp paired-end reads are 
sampled from both ends of these selected fragments.

Please refer to the paper snyder_rnaseq_review.pdf in the folder for
detailed RNA sequencing process.


## <a name="compilation"></a> Compilation & Installation

To compile simulator, simply run
   
    make

To install, simply put the simulator directory in your environment's PATH
variable.

### Prerequisites

C++ is required to be installed. 


## <a name="usage"></a> Usage
