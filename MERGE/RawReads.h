#ifndef RAWREADS_H
#define RAWREADS_H
#include<set>
#include<list>
#include<vector>
#include<iostream>
#include<algorithm>
#include<string>
#include <sstream>
#include <float.h>
#include <math.h>
#include <assert.h>   
#include <fstream>
#include <ctype.h>
#include <climits>
#include "TOK.h"
#include "TextIO.h"
#include "GenomeSegment.h"

class RawReads{
	public:
	void sep_gtfs(std::string fnamestr, std::string charname);
	
	void sep_gtfs_by_chrs(std::string fnamestr,std::string out_dir);
	void sep_reads(std::string fnamestr, std::string charname);
	void sep_reads_by_chrs(std::string out_dir , std::string fnamestr);
	void paired_reads_by_chrs(std::string out_dir );	
	void sep_reads_by_genes(std::string out_dir , std::string gtfname, std::string bedname);
	void bed2Fasta(std::string fname, std::vector<std::string> chr_dir);	
	void simulated_bed2Sam(std::string infname, std::string outfname);	
	void bed2Fasta2(std::string fname, std::vector<std::string> chr_dir);			
	void bed2Fasta_paired(std::string fname, std::vector<std::string> chr_dir);	
	void pre_process(std::string bedfname, std::string out_dir);
			
	std::string reverse_string(std::string in_str);
};

#endif