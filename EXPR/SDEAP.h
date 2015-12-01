#ifndef SDEAP_H
#define SDEAP_H
#include<map>
#include<set>
#include<list>
#include<vector>
#include<iostream>
#include<string>
#include <sstream>
#include <float.h>
#include <math.h>
#include "Merge.h"
#include "TextIO.h"
#include "TOK.h"
#include "RawReads.h"
#include "SPGraph.h"
#include "GenomeSegment.h"
#include "ASM.h"
#include "gene.h"
#include "dirent.h" 
class SDEAP{
	public:
	
	void runSDEAP(std::vector<std::string> BedVec, std::string gname,std::string gtfdir ,std::string tmpdir , std::string outdir, std::vector< double> lib_sizevec, std::vector< int> lenvec);
	void runRC(std::vector<std::string> BedVec, std::string gname, std::string gtfdir, std::string tmpdir , std::string outdir, std::vector< double> lib_sizevec);
	void runCOV(std::vector<std::string> BedVec, std::string gname,std::string gtfdir,std::string tmpdir , std::string outdir,std::vector< int> lenvec );
	
	std::vector<double> calcSizeFactor(std::vector<std::string> fnamevec);
	std::vector<std::string> var(std::string, int sample_size, double min_mean, double max_cv2 , double min_cv2_rate);
	std::vector<std::string> toGeneFeatures(std::string predfname);
	std::vector<std::string> toPCAdata(std::string predfname);	
	std::vector<std::string> 	toGeneRC(std::string indir, std::vector<std::string> gnamevec);
	void readSDEAP_Res(std::string allfname, std::string asfname);
		
	std::vector<double> Normalize4Dexus(std::vector<std::string> fnamevec);
	void cat_all(std::string outdir);
	int validCOV(std::string covfname); 	
	void sum_RCs(std::string glistfname,std::string outdir);
};

#endif