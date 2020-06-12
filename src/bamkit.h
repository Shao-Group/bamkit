#ifndef __BAMKIT_H__
#define __BAMKIT_H__

#include "hit.h"
#include <set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <iterator>

using namespace std;

typedef pair<pair<int, string>, pair<int, string> > pairPosCigar;
typedef pair<string, pairPosCigar> rcdIdentifier;

class bamkit
{
public:
	bamkit(const string &bamfile);
	~bamkit();

private:
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;

	int qcnt;			// single reads
	double qlen;		// single reads
	vector<int> ivec;	// insert size

    //evaluate aligners
    map<pair<string,int32_t>, bool > alignEvalMap;
    map<pair<string, int32_t>, pairPosCigar> hitIndexMap;
    map<rcdIdentifier, pair<string, int32_t> > hitIndexRevMap;
    set<rcdIdentifier> alignPairSet;

public:
	int solve_count();
	int solve_strand();
	int solve_fragment();;
	int ts2XS(const string &file);
	int name2to1(const string &file);
    int alignPairEval(const string &groundtruth);
    int bridgeEval(const string &alignerBam, const string &groundTruthBam, const string &annotation);
    int addXS(const string &file);
    int splitByEnd(const string &file1, const string &file2);
    int filter2ndAlign(const string &file);
    int splitSinglePaired(const string &file1, const string &file2);

private:
    int alignedPairs();
};

#endif
