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
