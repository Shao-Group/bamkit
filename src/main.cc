/*
Part of bamkit 
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "config.h"
#include "bamkit.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));

	if(argc != 3 && argc != 4)
	{
		printf("usage: \n");
		printf(" %s count <bam-file>\n", argv[0]);
		printf(" %s strand <bam-file>\n", argv[0]);
		printf(" %s fragment <bam-file>\n", argv[0]);
		printf(" %s ts2XS <in-bam-file> <out-bam-file>\n", argv[0]);
		printf(" %s name2to1 <in-bam-file> <out-bam-file>\n", argv[0]);
		return 0;
	}

	if(string(argv[1]) == "count")
	{
		bamkit bk(argv[2]);
		bk.solve_count();
	}

	if(string(argv[1]) == "strand")
	{
		bamkit bk(argv[2]);
		bk.solve_strand();
	}

	if(string(argv[1]) == "fragment")
	{
		bamkit bk(argv[2]);
		bk.solve_fragment();
	}

	if(string(argv[1]) == "ts2XS")
	{
		bamkit bk(argv[2]);
		bk.ts2XS(argv[3]);
	}

	if(string(argv[1]) == "name2to1")
	{
		bamkit bk(argv[2]);
		bk.name2to1(argv[3]);
	}

	return 0;
}
