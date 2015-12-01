#ifndef SIM_H
#define SIM_H
#include <map>
#include <set>
#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <float.h>
#include <math.h>
#include <limits>
#include "TOK.h"
#include "TextIO.h"
#include "gene.h"

class Sim{
	public:
	std::vector<std::string> getAllGeneList(std::string exprfname, std::string gtfdir, std::string idfname);
	std::map< std::string , double> calc_RC(std::string exprfname);	
	std::map< std::string , double> calc_RC_sp();
	
	std::vector<std::string> binary_up_RC(std::string glistfname,std::string spfname);
	std::vector<std::string> binary_down_RC(std::string glistfname,std::string spfname);		
	std::vector<std::string> binary_ds_RC(std::string glistfname,std::string spfname);
	std::vector<std::string> binary_EE_RC(std::string glistfname);	
	std::vector<std::string> binary_ALL_RC(std::string glistfname, std::vector<int> GroupVec, double dispersion);	
		std::vector<std::string> binary_ALL_RC(std::string glistfname, std::vector<int> GroupVec, double dispersion,int IDX);	
		
	std::map< std::string , std::string > read_transcripts_map(std::string gtfname , std::string gtfpath );	
	std::vector<std::string> createDSgenes(std::string exprfname, std::vector<std::string> testgeneVec, std::string idfname);	
	std::vector<std::string> simReads(std::map<std::string, std::string> ExprMap, std::map<std::string, std::string> BedMap, std::string tid, int read_num , int len, int bias);	
	std::vector<std::string> simReads_para(std::map<std::string, std::string> ExprMap, std::map<std::string, std::string> BedMap, std::string tid, int read_num , int len, int bias, int IDX);
	void simSamples(std::vector<int>  GroupVec ,	double lib_size ,double dispersion ,	std::string reads_outpath);
	void simSamples_para(std::vector<int>  GroupVec ,std::string rcfname,	double lib_size ,double dispersion ,	std::string reads_outpath, int IDX);
	void simSamples_test(std::vector<int>  GroupVec ,	double lib_size ,double dispersion ,	std::string reads_outpath );
	std::vector<std::string> shuffleExpr(std::vector<std::string> ExprVec);
	void baseMTX(std::string lenfname);
		
	std::vector<std::string> simReads(std::string ExprStr, std::string BedStr ,std::string tid, int read_num , int read_len, int bias);	
	void naive_2(int size_a, int size_b);			
	void naive_1(int size_a, int size_b);
	
	void update_RC();
	
	private:
			
		
};

#endif