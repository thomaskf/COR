#include "readUpdate.h"

//=====================================================
// For the class Nucl2Int
//=====================================================

// constructor
Nucl2Int::Nucl2Int() {
	// A->1, C->2, G->3, T->4
	memset(c, 0, 256);
	c[(int)'A'] = c[(int)'a'] = 1;
	c[(int)'C'] = c[(int)'c'] = 2;
	c[(int)'G'] = c[(int)'g'] = 3;
	c[(int)'T'] = c[(int)'t'] = 4;
}

//=====================================================
// For the class CoverageStat
//=====================================================

// constructor
CoverageStat::CoverageStat() {
	startMapPos = -1;
	consensusUpTo = 0;
	consensus = "";
}

void CoverageStat::addNucl(int nucl, int mpos) {
	int increaseNum = 10; // each time increase for 10 positions
	if (cov.size() == 0) {
		startMapPos = mpos;
	}
	int p = mpos - startMapPos;
	int p_size = p * 5;
	if (sum.size() <= p) {
		cov.resize(p_size + increaseNum*5, 0);
		sum.resize(p + increaseNum, 0);
	}
	cov[p_size + nucl]++;
	sum[p]++;
}

// compute consensus up to, but excluding, the map-pos
void CoverageStat::computeConsensus(int upto_mpos) {
	int i,j;
	int p;
	int max_cov;
	int max_nucl;
	int upToPos;
	char int2nucl[] = {'N','A','C','G','T'};

	if (startMapPos == -1)
		return;
	
	if (upto_mpos == -1) {
		upToPos = sum.size();
	} else {
		upToPos = upto_mpos - startMapPos;
	}
	if (consensusUpTo < upToPos) {
		for (i=consensusUpTo; i<upToPos; i++) {
			if (i < sum.size() && sum[i] >= COVER_THRES) {
				p = i*5;
				max_cov = cov[p];
				max_nucl = 0;
				for (j=1; j<=4; j++) {
					if (cov[p+j] > max_cov) {
						max_cov = cov[p+j];
						max_nucl = j;
					}
				}
				consensus.append(1,int2nucl[max_nucl]);
			} else {
				consensus.append(1,int2nucl[0]);
			}
		}
		consensusUpTo = upToPos;
	}
}

int CoverageStat::consensusUptoMapPos() {
	return consensusUpTo + startMapPos;
}

void CoverageStat::show(int frPos) {
	int i,j;
	if (startMapPos == -1)
		return;
		
	cout << "Coverage:" << endl;
	for (i = 0; i < cov.size(); i+=5) {
		int p = startMapPos + i/5;
		if (frPos == -1 || startMapPos >= frPos) {
			cout << p << " :";
			for (j = 0; j < 5; j++) {
				cout << " " << cov[i+j];
			}
			cout << " [" << sum[i/5] << "]" << endl;
		}
	}
	
	cout << "Consensus:" << endl;
	if (frPos == -1) {
		cout << consensus << endl;
	} else {
		int p = frPos - startMapPos;
		int len = consensus.length() - p;
		cout << consensus.substr(p, len) << endl;
	}
}

//=====================================================
// For the class Read
//=====================================================

// constructor
Read::Read(string& rd_name, string& rd_seq, string& rd_qualseq, vector<int>& mPoses, int endmpos, string& samseq) {
	name = rd_name;
	seq = rd_seq;
	qualseq = rd_qualseq;
	mapPoses.clear();
	mapPoses.insert(mapPoses.begin(), mPoses.begin(), mPoses.end());
	endMapPos = endmpos;
	SAMseq = samseq;
	startMapPos = -1;
	for (int i=0; i<mapPoses.size(); i++) {
		if (mapPoses[i] != -1) {
			startMapPos = mapPoses[i];
			break;
		}
	}
	toDiscard = false;
}

Read::Read(string& samseq) {
	SAMseq = samseq;
	toDiscard = true;
}

// correct this read according to the consensus
// return true if the read is covered by the consensus sequence
bool Read::correctFromConsensus(CoverageStat& coverageStat, bool& changed) {
	int i;
	int mpos, relpos;
	changed = false;

	if (toDiscard)
		return true;
		
	if (endMapPos < coverageStat.consensusUptoMapPos()) {
		for (i=0; i<mapPoses.size(); i++) {
			mpos = mapPoses[i];
			if (mpos != -1) {
				relpos = mpos - coverageStat.startMapPos;
				if (coverageStat.consensus[relpos] != 'N') {
					if (seq[i] != coverageStat.consensus[relpos]) {
						seq[i] = coverageStat.consensus[relpos];
						qualseq[i] = QUAL_VALUE_CHANGED_BASE + '!';
						changed = true;
					}
				}
			}
		}
		return true;
	}
	return false;
}


// output the read to the SAM file
void Read::outputSAM(ofstream& fout) {
	if (toDiscard) {
		fout << SAMseq << endl;
		return;
	}
	vector<string> token;
	tokenizer(SAMseq, "\t", &token);
	for (int i=0; i<token.size(); i++) {
		if (i>0)
			fout << "\t";
		if (i==9) {
			// output the updated sequence
			fout << seq;
		} else if (i==10) {
			// output the updated quality sequence
			fout << qualseq;
		} else {
			fout << token[i];
		}
	}
	fout << endl;
}

// check the seq according to the cigar string
// get the mapping positions of each character along the sequence
// return false if the cigarStr is not valid
bool checkSeq(string& seq, string& cigarStr, int readId, int mapPos, int& frontTrim, int& endTrim, string& qualseq, vector<int>& mPoses, int& endPos) {
	int i,j,k,l;
	int seqI;
	
	i=0;
	seqI = 0;
	frontTrim = 0;
	endTrim = 0;
	mPoses.clear();
	for (j=0; j<cigarStr.length(); j++) {
		// should begin with a number, and followed by a char
		if (isdigit(cigarStr[j]))
			continue;
		else {
			if (j==i) {
				// cout << "Warning! Invalid cigar string: " << cigarStr << endl;
				return false;
			}
			k = atoi(cigarStr.substr(i,j-i).c_str());
			switch (cigarStr[j]) {
				case 'I':
					// inserts.push_back(new Insert(mapPos+seqI,seqI,seq.substr(seqI,k),readId));
				case 'S':
					// seq.erase(seqI,k);
					// qualseq.erase(seqI,k);
					if (seqI == 0)
						frontTrim = k;
					if (j==cigarStr.length()-1)
						endTrim = k;
					for (l=0; l<k; l++)
						mPoses.push_back(-1);
					break;
				case 'D':
				case 'N':
					// seq.insert(seqI,k,'-');
					// qualseq.insert(seqI,k,0);
					seqI += k;
					break;
				case 'M':
					for (l=0; l<k; l++)
						mPoses.push_back(mapPos + seqI + l);
					endPos = mapPos + seqI + k - 1;
					seqI += k;
					break;
				case 'P':
				case 'H':
					// do nothing
					break;
			}
			i=j+1;
		}
	}
	return true;
}

// update the coverage statistics
void updateStat(string& seq, string& qualseq, vector<int>& mPoses, CoverageStat& coverageStat, Nucl2Int& nucl2int) {

	int qual;
	int nucl;
	int i,startpos;
	int currRelPos = 0;
	
	// find the first position i where the mPoses[i] is not -1
	startpos = 0;
	while (startpos<mPoses.size() && mPoses[startpos] == -1) {
		startpos++;
	}
	
	if (startpos >= mPoses.size()) {
		// all mPoses[i] are -1
		return;
	}
	
	for (i=startpos; i<qualseq.size(); i++) {
		if (mPoses[i] != -1) {
			qual = qualseq[i] - '!';
			if (qual >= QV_THRES) {
				// the base's quality is good enough
				nucl = (int) nucl2int.c[seq[i]];
				// add the nucleotide information to the coveragestat
				coverageStat.addNucl(nucl, mPoses[i]);
			}
		}
	}
}

// update the reads according to the consensus
// then output the reads to a fastq file
void updateReads(char* samFile, char* outFQFile) {
	string aline;
	int mpos;
	int mapq;
	int flag;
	int indx;
	int frontTrim;
	int endTrim;
	int endPos;
	int endPosUnreportRead = -1;
	Nucl2Int nucl2int;
	string cigar;
	string seq;
	string qualseq;
	vector<string> token;
	SamBamFileHander fileHander;
	string preSeqName = "";
	string currSeqName;
	vector<int> mPoses;
	list<Read> reads;
	CoverageStat coverageStat;
	ofstream fout;
	bool toDiscard;
	bool changed;
	indx = 0;
	list<Read>::iterator itr_rd;
	int totNumReads = 0;
	int changedReads = 0;
	int ignoredReads = 0;
	
	// open the output file
	fout.open(outFQFile);
	
	// store the unprocessed reads
	fileHander.openFile(samFile);
	while (fileHander.getNextSeq(aline)) {
		indx++;
		if (aline.length() > 0 && aline[0]!='@') {
			tokenizer(aline, "\t", &token);
			if (token.size() > 9) {
				currSeqName = token[0];
				toDiscard = false;
				totNumReads++;
				flag = atoi(token[1].c_str());
				if ((flag&2048)!=0) {
					// this is a supplementry alignment
					// discard it
					toDiscard = true;
				}
				mpos = atoi(token[3].c_str());
				if (mpos==0) {
					// this read is not aligned
					toDiscard = true;
				}
				mpos = mpos-1; // change to zero-based
				mapq = atoi(token[4].c_str());
				cigar = token[5];
				if (cigar == "*") {
					// this read is not aligned
					toDiscard = true;
				}
				seq = token[9];
				if (mapq < MAPQ_THRES) {
					// the mapping quality is too low, thus this read is ignored
					// cout << aline << endl;
					toDiscard = true;
				}
				qualseq = token[10];
				if (!toDiscard && checkSeq(seq, cigar, indx, mpos, frontTrim, endTrim, qualseq, mPoses, endPos)) {
				
					if ((double)(frontTrim + endTrim) / seq.length() > MAP_TRIMLEN_RATIO_THRES) {
						// the trimmed region is too long, thus this read is ignored
						toDiscard = true;
					}

					if (!toDiscard) {
						reads.push_back(Read(currSeqName, seq, qualseq, mPoses, endPos, aline));
					
						if (endPosUnreportRead == -1)
							endPosUnreportRead = endPos;

						// compute consensus up to, but excluding, the map-pos
						coverageStat.computeConsensus(mpos);
					
						// output the corrected reads if there exists
						int numReadToRemove = 0;
						for (itr_rd = reads.begin(); itr_rd != reads.end(); itr_rd++) {
							if (itr_rd->correctFromConsensus(coverageStat, changed)) {
								itr_rd->outputSAM(fout);
								numReadToRemove++;
								if (changed)
									changedReads++;
							} else {
								break;
							}
						}
					
						// remove the reads which have been corrected and reported to the FASTQ file
						for (int i=0; i<numReadToRemove; i++)
							reads.pop_front();

						// update the statistics
						updateStat(seq, qualseq, mPoses, coverageStat, nucl2int);
					}
					
				} 
				
				if (toDiscard) {
					ignoredReads++;
					reads.push_back(Read(aline));
				}
			}
		} else {
			// it is comment
			fout << aline << endl;
		}
	}
	fileHander.closeFile();
	
	// compute consensus to the end
	coverageStat.computeConsensus(-1);

	// output the rest of the reads
	for (itr_rd = reads.begin(); itr_rd != reads.end(); itr_rd++) {
		if (itr_rd->correctFromConsensus(coverageStat, changed)) {
			itr_rd->outputSAM(fout);
			if (changed)
				changedReads++;
		} else {
			cout << "Error! This read " << itr_rd->name << " cannot be corrected!!!" << endl;
		}
	}
	
	// close the output file
	fout.close();
	
	// print out the statistics
	cout << "Total number of reads: " << totNumReads << endl;
	cout << "Total number of reads ignored: " << ignoredReads << endl;
	cout << "Total number of reads updated: " << changedReads << endl;
}

