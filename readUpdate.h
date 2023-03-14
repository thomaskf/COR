#ifndef _READ_UPDATE_
#define _READ_UPDATE_

#include<string>
#include<list>
#include<fstream>
#include<iostream>
#include<cstdlib>
#include "fileHandler.h"
#include "mylib.h"
#include "definitions.h"

using namespace std;

class Nucl2Int {
public:
	char c[256];
	// constructor
	Nucl2Int();
};

class CoverageStat {
public:
	int startMapPos;
	int consensusUpTo; // the consensus sequence has been computed from the beginning up to, but excluding, this pos
	vector<int> cov; // coverage for each character along the positions
	vector<int> sum; // total coverage  along the positions
	string consensus; 	// consensus sequence
	// constructor
	CoverageStat();
	// add a read with the corresponding nucleotide to the position
	void addNucl(int nucl, int mpos);
	// compute consensus up to, but excluding, the pos
	void computeConsensus(int uptoPos);
	// output the mapPos position up to
	int consensusUptoMapPos();
	// show
	void show(int frPos = -1);
};

class Read {
public:
	string name;
	string seq;
	string qualseq;
	vector<int> mapPoses;
	int startMapPos;
	int endMapPos;
	bool toDiscard;
	
	string SAMseq;

	// constructor
    Read(string& rd_name, string& rd_seq, string& rd_qualseq, vector<int>& mPoses, int endmpos, string& samseq);
    Read(string& samseq);
    
    // correct this read according to the consensus
    // return true if the read is covered by the consensus sequence
    bool correctFromConsensus(CoverageStat& coverageStat, bool& changed);

	// output the read to the SAM file
	void outputSAM(ofstream& fout);
};

// update the reads according to the consensus
void updateReads(char* samFile, char* outFQFile);

#endif