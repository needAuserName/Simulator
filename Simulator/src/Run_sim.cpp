#include <fstream>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <queue>
#include <stack>
#include <sstream>

using namespace std;

// Transform long to string
string itostr(long t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}

int Round = 5;
int CopyNumCnt = 3;
int CopyNum[3] = {1,2,5};

int main(int argc, char* argv[])
{
	string comd, resultPath, outpath;
	int i;

	// 2. parse frag into sam
	resultPath = "/data/homes/yan/Transcript_Assembly/Simulation/simulation_gap/src/single_gene_example/Gap_investigation/Long_transcript/";
	string filename1, filename2, filepath1;
// 	for (int copyIndex = 0; copyIndex < CopyNumCnt; copyIndex++)
// 	{
// 		for (int round = 1; round <= Round; round++)
// 		{
// 			filename1 = resultPath + "copyNum=" + itostr(CopyNum[copyIndex]) + "/round" + itostr(round) + "/allfraginfo.txt";
// 			filename2 = resultPath + "test.gtf";
// 			filepath1 = resultPath + "copyNum=" + itostr(CopyNum[copyIndex]) + "/round" + itostr(round) + "/allfrag.sam";
// 			comd = "/data/homes/yan/Transcript_Assembly/Simulation/simulation_gap/src/single_gene_example/src/GenerateSAMGTF " + filename1 + " " + filename2 + " " + filepath1;
// 			i = system(comd.c_str());
// 
// 			filename1 = resultPath + "copyNum=" + itostr(CopyNum[copyIndex]) + "/round" + itostr(round) + "/fraginfo_selected.txt";
// 			filename2 = resultPath + "test.gtf";
// 			filepath1 = resultPath + "copyNum=" + itostr(CopyNum[copyIndex]) + "/round" + itostr(round) + "/selectedfrag.sam";
// 			comd = "/data/homes/yan/Transcript_Assembly/Simulation/simulation_gap/src/single_gene_example/src/GenerateSAMGTF " + filename1 + " " + filename2 + " " + filepath1;
// 			i = system(comd.c_str());
// 		}
// 	}

	// 3. parse fragments
	for (int copyIndex = 0; copyIndex < CopyNumCnt; copyIndex++)
	{
		for (int round = 1; round <= Round; round++)
		{
			comd = "mkdir " + resultPath + "copyNum=" + itostr(CopyNum[copyIndex]) + "/round" + itostr(round) + "/selectedfrag";
			i = system(comd.c_str());
			comd = "./parsefrag " + resultPath + "copyNum=" + itostr(CopyNum[copyIndex]) + "/round" + itostr(round) + "/selectedfrag.sam 1 noDB " + resultPath + "copyNum=" + itostr(CopyNum[copyIndex]) + "/round" + itostr(round) + "/selectedfrag/";
			i = system(comd.c_str());
		}
	}

	// 4. assembly (feed annotation)


	return 0;
}