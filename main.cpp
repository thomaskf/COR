#include "readUpdate.h"
#include "mylib.h"
#include <math.h>
#include <time.h>
#include <fstream>
#include <vector>

#define VERSION "c1.0"

void showHelpMenu(char** argv) {
    cerr << "COR (Consensus Overrides Reads) aims to reduce the sequencing errors and remove the heterozygous information on the reads in order to reduce the complexity in the subsequent assembly process. COR updates the read sequences, which should be from the same species, according to the consensus of the read alignments against a reference sequence." << endl;
    cerr << endl;
    cerr << "Syntax:" << endl;
    cerr << "  " << argv[0] << " [in: sorted SAM file] [out: sorted SAM file]" << endl;
    cerr << endl;
    cerr << "Please note that the program does not work if the input SAM file is not sorted." << endl;
    cerr << endl;
}

int main(int argc, char** argv) {
    
	// display version number
	cout << "Version " << VERSION << endl;
    // cout << "This version considers substrings with two snp positions" << endl;
	cout << endl;
	
    // check the syntax
    if (argc < 3) {
        showHelpMenu(argv);
        exit(1);
    }

    char* inputFile = argv[1];
    char* outputFile = argv[2];
    int i,j;

	// update the reads according to the consensus
	updateReads(inputFile, outputFile);
}
