
#ifndef GENE_H
#define GENE_H
#include<map>
#include<set>
#include<list>
#include<vector>
#include<iostream>
#include<string>
#include <sstream>
#include <float.h>
#include <math.h>
#include <limits.h>
#include "TOK.h"
#include "GenomeSegment.h"
#include "SPGraph.h"
#include "ASM.h"
#include "RawReads.h"
class gene{
	public:
  GenomeSegment gs;
  SPGraph spg ; 
  SPGraph spg_rc ; 
  
  ASM D;
  SPGraph build_from_GTF(std::string gtf_gname);
  std::vector<std::string> init(std::string gtf_gname,std::string bed_fname,std::string gname, double scale);	
  std::vector<std::string> init_SPG(std::string gtf_gname,std::string bed_fname,std::string gname,std::string tmpdir,int reads_len, double scale);	
	
	std::vector<std::string> init_Hybrid(std::string gtf_gname,std::string bed_fname,std::string gname,std::string tmpdir,int reads_len, double scale);
	
	std::vector<std::string> init_Hybrid_RC(std::string gtf_gname,std::string bed_fname,std::string gname,std::string tmpdir,int reads_len, double scale);
	void init_SPG_exon(std::string gtf_gname,std::string bed_fname,std::string gname,std::string tmpdir,int reads_len, double scale); 
  int  multiGene(std::string gtf_gname);
  void create_coverage(std::string gtf_gname, std::string bed_fname);	
  void map_one_frag_bed(std::string read1_str , std::string read2_str);
  void map_one_read_bed_toSPG(std::string read1_str, double scale);	
  void map_one_read_bed(std::string read1_str, double scale);	
  void map_one_bed_file(std::string fname, double scale);	
  void map_one_bed_file_toSPG(std::string fname, double scale);
	
	
};

#endif