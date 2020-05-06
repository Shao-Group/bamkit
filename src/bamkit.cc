#include <cstdio>
#include <cassert>
#include <sstream>

#include "config.h"
#include "bamkit.h"

bamkit::bamkit(const string &bamfile)
{
    sfn = sam_open(bamfile.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	qlen = 0;
	qcnt = 0;
}

bamkit::~bamkit()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
}

int bamkit::solve_count()
{
	int maxisize = 500;
	ivec.clear();
	ivec.assign(maxisize, 0);
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > MAX_NUM_CIGAR) continue;									// ignore hits with more than 7 cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t, 1);

		qlen += ht.qlen;
		qcnt += 1;

		if(ht.isize <= 0 || ht.isize >= maxisize) continue;

		ivec[ht.isize]++;
	}

	int icnt = 0;
	double iave = 0;
	double idev = 0;
	for(int i = 1; i < ivec.size(); i++)
	{
		icnt += ivec[i];
		iave += ivec[i] * i;
	}
	iave = iave / icnt;

	for(int i = 1; i < ivec.size(); i++)
	{
		idev += (i - iave) * (i - iave) * ivec[i];
	}
	idev = sqrt(idev / icnt);

	printf("aligned reads = %d aligned base pair = %.0lf average read length = %.2lf insert size = %.2lf +- %.2lf\n", qcnt, qlen, qlen / qcnt, iave, idev);

	return 0;
}

int bamkit::solve_strand()
{
	int first = 0;
	int second = 0;
	int cnt = 0;
	int n = 100000;

	library_type = FR_FIRST;
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;			// read is not mapped
		if((p.flag & 0x8) >= 1) continue;			// mate is note mapped
		if((p.flag & 0x100) >= 1) continue;			// secondary alignment
		if(p.n_cigar > MAX_NUM_CIGAR) continue;

		hit ht(b1t, 6);
        //cout << ht.hi << endl;
		if(ht.xs == '.') continue;

		if(ht.strand == '+' && ht.xs == '+') first++;
		if(ht.strand == '-' && ht.xs == '-') first++;
		if(ht.strand == '+' && ht.xs == '-') second++;
		if(ht.strand == '-' && ht.xs == '+') second++;
		
		//printf("xs = %c, strand = %c\n", ht.xs, ht.strand);

		cnt++;
		if(cnt >= n) break;
	}

	string type = "unstranded";
	if(cnt >= 0.8 * n && first >= 0.8 * cnt) type = "first";
	if(cnt >= 0.8 * n && second >= 0.8 * cnt) type = "second";

	printf("samples = %d, first = %d, second = %d, library = %s\n", cnt, first, second, type.c_str());
	return 0;
}

int bamkit::solve_fragment()
{
	int maxisize = 500;
	ivec.clear();
	ivec.assign(maxisize, 0);
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > MAX_NUM_CIGAR) continue;									// ignore hits with more than 7 cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar != 1) continue;

		hit ht(b1t, 1);

		qlen += ht.qlen;
		qcnt += 1;

		if(ht.isize <= 0 || ht.isize >= maxisize) continue;

		ivec[ht.isize]++;
	}

	int icnt = 0;
	double iave = 0;
	double idev = 0;
	for(int i = 1; i < ivec.size(); i++)
	{
		icnt += ivec[i];
		iave += ivec[i] * i;
	}
	iave = iave / icnt;

	for(int i = 1; i < ivec.size(); i++)
	{
		idev += (i - iave) * (i - iave) * ivec[i];
	}
	idev = sqrt(idev / icnt);

	printf("aligned reads = %d aligned base pair = %.0lf average read length = %.2lf insert size = %.2lf +- %.2lf\n", qcnt, qlen, qlen / qcnt, iave, idev);

	return 0;
}

int bamkit::ts2XS(const string &file)
{
	samFile *fout = sam_open(file.c_str(), "w");

	int f = sam_hdr_write(fout, hdr);
	if(f < 0) printf("fail to write header to %s\n", file.c_str());
	if(f < 0) exit(0);

    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		uint8_t *p = bam_aux_get(b1t, "ts");

		if((p) && (*p) == 'A')
		{
			char XS = '.';
			char ts = bam_aux2A(p);
			if(ts == '+' && (((b1t)->core.flag) & 0x10) <= 0) XS = '+';
			if(ts == '+' && (((b1t)->core.flag) & 0x10) >= 1) XS = '-';
			if(ts == '-' && (((b1t)->core.flag) & 0x10) <= 0) XS = '-';
			if(ts == '-' && (((b1t)->core.flag) & 0x10) >= 1) XS = '+';

			f = bam_aux_append(b1t, "XS", 'A', sizeof(XS), (uint8_t *) &XS);
			if(f != 0) printf("fail to append XS\n");
			if(f != 0) exit(0);
		}

		f = sam_write1(fout, hdr, b1t);
		if(f < 0) printf("fail write alignment to %s\n", file.c_str());
		if(f < 0) exit(0);
	}

	sam_close(fout);
	return 0;
}

int bamkit::name2to1(const string &file)
{
	BGZF *fout = bgzf_open(file.c_str(), "w");
	int f = bam_hdr_write(fout, hdr);
	if(f < 0) printf("fail to write header to %s\n", file.c_str());
	if(f < 0) exit(0);

    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		int l = b1t->core.l_qname - b1t->core.l_extranul - 1;
		char *qname = bam_get_qname(b1t);
		assert(l >= 2);
		assert(qname[l - 2] == '.');
		if(qname[l - 1] == '2') qname[l - 1] = '1';
		f = bam_write1(fout, b1t);
		if(f < 0) printf("fail write alignment to %s\n", file.c_str());
		if(f < 0) exit(0);
	}

	bgzf_close(fout);
	return 0;
}

int bamkit::alignPairEval(const string &groundtruth)
{
    bamkit gt(groundtruth);

    alignedPairs();
    gt.alignedPairs();
    
    set<rcdIdentifier> commonSet, unalignedSet, wrongSet;
    set_intersection(gt.alignPairSet.begin(), gt.alignPairSet.end(), alignPairSet.begin(), alignPairSet.end(), inserter(commonSet, commonSet.begin()));
    set_difference(gt.alignPairSet.begin(), gt.alignPairSet.end(), alignPairSet.begin(), alignPairSet.end(), inserter(unalignedSet, unalignedSet.begin()));
    set_difference(alignPairSet.begin(), alignPairSet.end(), gt.alignPairSet.begin(), gt.alignPairSet.end(), inserter(wrongSet, wrongSet.begin()));


    FILE * commonFile = fopen("common.txt", "w");
    FILE * wrongFile = fopen("wrong.txt", "w");
    for(auto it = commonSet.begin(); it != commonSet.end(); it++)
    {
        alignEvalMap[hitIndexRevMap[*it]] = true;
        fprintf(commonFile, "%s HI:%d\n", hitIndexRevMap[*it].first.c_str(), hitIndexRevMap[*it].second);
    }

    for(auto it = wrongSet.begin(); it != wrongSet.end(); it++)
    {
        fprintf(wrongFile, "%s HI:%d\n", hitIndexRevMap[*it].first.c_str(), hitIndexRevMap[*it].second);
    }
    fclose(commonFile);
    fclose(wrongFile);
    
    uint32_t common = commonSet.size(),unaligned = unalignedSet.size(), wrong = wrongSet.size();
    uint32_t totalGT = gt.alignPairSet.size(), totalAligner = alignPairSet.size();    
    float sensitivity = 1.0*common/totalGT;
    float precision = 1.0*common/totalAligner;
    printf("Total of ground truth: %d, Total of aligner: %d, Unaligned: %d, True alignment: %d, False alignment:%d\n",totalGT, totalAligner, unaligned, common, wrong);
    printf("Sensitivity: %.4f, Precision: %.4f\n", sensitivity, precision);
    return 0;
}

int bamkit::alignedPairs()
{
    hitIndexMap.clear();
    hitIndexRevMap.clear();
    alignPairSet.clear();
    library_type = FR_SECOND;//flux simulated data is FRsecond
    while(sam_read1(sfn,hdr,b1t) >= 0)
    {
        string qname = bam_get_qname(b1t);
        uint32_t* cigar = bam_get_cigar(b1t);
        bam1_core_t &p = b1t->core;
        
        if((p.flag & 0x4) >= 1) continue;
        //if((p.flag & 0x100) >= 1) continue;

        string cigarStr = "";
        for(int i=0; i < p.n_cigar;++i)
        {
            int icigar = cigar[i];
            cigarStr=cigarStr + to_string(bam_cigar_oplen(icigar)) + bam_cigar_opchr(icigar);
        }
        
        hit ht(b1t, 4);
        pair<string, uint32_t> key = make_pair(qname, ht.hi);
        alignEvalMap[key] = false;
        if(((p.flag & 0x40) >= 1 && ht.strand == '+') || ((p.flag & 0x80) >= 1 && ht.strand == '-'))//first segment
        {
            if(hitIndexMap.find(key) == hitIndexMap.end())
            {
                hitIndexMap[key] = make_pair(make_pair(p.pos, cigarStr), make_pair(-1, ""));
            }
            else
            {
                hitIndexMap[key].first = make_pair(p.pos, cigarStr);
            }
        }
        else if(((p.flag & 0x40) >= 1 && ht.strand == '-') || ((p.flag & 0x80) >= 1 && ht.strand == '+'))
        {
            if(hitIndexMap.find(key) == hitIndexMap.end())
            {
                hitIndexMap[key] = make_pair(make_pair(-1, ""), make_pair(p.pos, cigarStr));
            }
            else
            {
                hitIndexMap[key].second = make_pair(p.pos, cigarStr);
            }
        }
    }

    for(auto it = hitIndexMap.begin(); it != hitIndexMap.end(); it++)
    {
        //printf("%s HI:%d [(%d, %s), (%d, %s)]\n", (it->first).first.c_str(), (it->first).second, (it->second).first.first, (it->second).first.second.c_str(), (it->second).second.first, (it->second).second.second.c_str());
        string qname = (it->first).first;
        hitIndexRevMap[make_pair(qname, it->second)] = it->first;
        alignPairSet.insert(make_pair(qname, it->second));
    }
}


int bamkit::bridgeEval(const string &alignerBam, const string &groundTruthBam, const string &annotation)
{
    bamkit aligner(alignerBam);
    aligner.alignPairEval(groundTruthBam);
    alignEvalMap = aligner.alignEvalMap;
    /*for(auto it = alignEvalMap.begin(); it != alignEvalMap.end(); it++)
    {
        printf("[%s, HI:%d]\n", (it->first).first.c_str(), (it->first).second);
    }
    cout << "lala" << endl;*/

    map< string, set< pair<uint32_t,uint32_t> > > exonMap;
    ifstream fin(annotation);

    string line;
    while(getline(fin, line))
    {
        stringstream liness(line);

        string chrid;
        string source;
        string type;
        getline(liness, chrid, '\t');
        getline(liness, source, '\t');
        getline(liness, type, '\t');
        if(type == "exon")
        {
            string startss, endss, score, strand, frame, info;
            getline(liness, startss, '\t');
            getline(liness, endss, '\t');
            getline(liness, score, '\t');
            getline(liness, strand, '\t');
            getline(liness, frame, '\t');
            getline(liness, info, '\t');
            
            uint32_t start = stoi(startss);
            uint32_t end = stoi(endss);
            
            stringstream infoss(info);
            string tid;
            while(getline(infoss,tid, ' ') && tid != "transcript_id");
            getline(infoss, tid, ' ');
            if(tid[0] == '\"')
            {
                tid.erase(tid.begin());
                if(tid.back() == ';')
                    tid.erase(tid.end()-2, tid.end());
            }

            if(exonMap.find(tid) == exonMap.end())
            {
                set< pair<uint32_t, uint32_t> > s;
                s.insert(make_pair(start-1, end));//0-based, half-open
                exonMap[tid] = s;
                //printf("%s: (%d, %d)\n", tid.c_str(), start, end);
            }
            else
            {
                exonMap[tid].insert(make_pair(start-1,end));//0-based, half-open
                //printf("%s: (%d, %d)\n", tid.c_str(), start, end);
            }
        }
    }


    //build ground truth bridge
    //assume the ground-truth bam is FR-first: R1+,R2-
    //fragment:[pos, rpos)
    bamkit gt(groundTruthBam);
    map< string, pair<uint32_t, string> > bridgeMap;//qname->(start, cigar)
    map< string, pair<uint32_t, uint32_t> > fragmentMap;//qname->(start,end)
    string qname;
    while(sam_read1(gt.sfn, gt.hdr, gt.b1t) >= 0)
    {
        qname = bam_get_qname(gt.b1t);
        bam1_core_t &p = gt.b1t->core;
        if((p.flag & 0x40) >= 1)
        {
            if(fragmentMap.find(qname) == fragmentMap.end())
                fragmentMap[qname] = make_pair(p.pos, 0);
            else
                fragmentMap[qname].first = p.pos;
        }
        else if((p.flag & 0x80) >= 1)
        {
            uint32_t frEnd = p.pos+(int32_t)bam_cigar2rlen(p.n_cigar, bam_get_cigar(gt.b1t));
            if(fragmentMap.find(qname) == fragmentMap.end())
                fragmentMap[qname] = make_pair(0,frEnd);
            else
                fragmentMap[qname].second = frEnd;
        }

        //printf("%s: (%d, %d)\n", qname.c_str(), fragmentMap[qname].first, fragmentMap[qname].second);
    }

    for(auto it = fragmentMap.begin(); it != fragmentMap.end(); it++)
    {   
        stringstream qnamess(it->first);
        string chr, locus, tr;
        getline(qnamess, chr, ':');
        getline(qnamess, locus, ':');
        getline(qnamess, tr, ':');
        //cout << tr << endl;
        uint32_t brStart = (it->second).first, brEnd = (it->second).second;
        uint32_t p1 = brStart, p2 = brStart;
        //printf("Paired read: (%d, %d)\n", brStart, brEnd);
        string newCigar = "";
        for(auto it = exonMap[tr].begin(); it != exonMap[tr].end() && p2<brEnd; it++)
        {
            //printf("Exon: (%d, %d)\n", it->first, it->second);
            if(brStart> it->second) continue;

            p1 = max(p2, it->first);
            if(p1-p2>0)
                newCigar = newCigar+ to_string(p1-p2) + "N";

            p2 = min(brEnd, it->second);
            if(p2-p1>0)
                newCigar = newCigar+ to_string(p2-p1)+"M";
            //printf("p1: %d; p2: %d; %s\n", p1, p2, newCigar.c_str());
        }
        bridgeMap[it->first] = make_pair(brStart, newCigar);
        //printf("%s: (%d, %s)\n", (it->first).c_str(), brStart, newCigar.c_str());
    }
    
    uint64_t cntTotalTruth = fragmentMap.size(), cntBridged = 0, cntBridgedCorrect = 0;
    uint64_t unbridged = 0, trueBrFalseAl = 0, trueBrTrueAl = 0, falseBrFalseAl = 0, falseBrTrueAl = 0;
    uint64_t unBrTrueAl = 0, unBrFalseAl = 0;
    uint64_t trueBrUnal = 0, falseBrUnal = 0, unBrUnal = 0;
    uint32_t *cigar;
    FILE * matchFile = fopen ("matchFile.txt","w");
    FILE * mismatchFile = fopen("mismatchFile.txt", "w");
    while(sam_read1(sfn, hdr, b1t) >= 0)
    {
        qname = bam_get_qname(b1t);
        bam1_core_t &p = b1t->core;
        hit ht(b1t, 4);
        pair<string, int32_t> key = make_pair(qname, ht.hi);
        //printf("[%s, HI:%d]\n", key.first.c_str(), key.second);
        
        if((p.flag & 0x4) >= 1) continue;
        //if((p.flag & 0x100) >= 1) continue;
        
        if(p.mpos != 0) 
        {
            if(((p.flag & 0x40) >= 1 && ht.strand == '+') || ((p.flag & 0x80) >= 1 && ht.strand == '-'))//first segment
            {
                if(alignEvalMap.find(key) == alignEvalMap.end())
                    unBrUnal++;
                else if(alignEvalMap[key])
                    unBrTrueAl++;
                else
                    unBrFalseAl++;
                unbridged++;
            }
            continue;//not bridged
        }
        cntBridged++;

        //cout << p.n_cigar << endl;
        cigar = bam_get_cigar(b1t);
        string cigarStr = "";
        for(int i=0; i < p.n_cigar;++i)
        {
            int icigar = cigar[i];
            cigarStr=cigarStr + to_string(bam_cigar_oplen(icigar)) + bam_cigar_opchr(icigar);
        }

        if(bridgeMap[qname].first == p.pos && bridgeMap[qname].second == cigarStr)
        {
            if(alignEvalMap.find(key) == alignEvalMap.end())
                trueBrUnal++;
            else if(alignEvalMap[key])
                trueBrTrueAl++;
            else
            {
                trueBrFalseAl++;
                fprintf(matchFile, "%s:\nGT:(%d, %s)\tCoral:(HI:%d, %d, %s)\n",qname.c_str(), bridgeMap[qname].first, bridgeMap[qname].second.c_str(), ht.hi, p.pos, cigarStr.c_str());
            }

            cntBridgedCorrect++;
            //fprintf(matchFile, "%s:\nGT:(%d, %s)\tCoral:(HI:%d, %d, %s)\n",qname.c_str(), bridgeMap[qname].first, bridgeMap[qname].second.c_str(), ht.hi, p.pos, cigarStr.c_str());
        }
        else
        {
            if(alignEvalMap.find(key) == alignEvalMap.end())
                falseBrUnal++;
            else if(alignEvalMap[key])
                falseBrTrueAl++;
            else
                falseBrFalseAl++;

            //fprintf(mismatchFile, "%s:\nGT:(%d, %s)\tCoral:(HI:%d, %d, %s)\n",qname.c_str(), bridgeMap[qname].first, bridgeMap[qname].second.c_str(), ht.hi, p.pos, cigarStr.c_str());
        }

    }
    fclose(matchFile);
    fclose(mismatchFile);

    float sensitivity = 1.0*cntBridgedCorrect/cntTotalTruth;
    float precision = 1.0*cntBridgedCorrect/cntBridged;
    printf("#pairs_in_ground_truth: %ld, #bridged_pairs: %ld, #correct_bridged_pairs:%ld\n", cntTotalTruth, cntBridged, cntBridgedCorrect);
    printf("#Unbridged_unalign: %ld, #True_bridged_unalign: %ld, #False_bridged_unalign: %ld\n", unBrUnal, trueBrUnal, falseBrUnal);
    printf("#Unbridged_true_align: %ld, #Unbridged_false_align: %ld\n", unBrTrueAl, unBrFalseAl);
    printf("#True_bridged_false_align: %ld, #True_bridged_true_align: %ld\n#False_bridged_false_align: %ld, #False_bridged_true_align: %ld\n", trueBrFalseAl, trueBrTrueAl, falseBrFalseAl, falseBrTrueAl);
    printf("Sensitivity: %.4f, Precision: %.4f\n",sensitivity, precision);
    return 0;
}

