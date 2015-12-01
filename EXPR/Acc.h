#ifndef ACC_H
#define ACC_H
#include<map>
#include<set>
#include<list>
#include<vector>
#include<iostream>
#include<string>
#include <sstream>
#include <float.h>
#include <math.h>
#include "TOK.h"
#include "TextIO.h"
#include "dirent.h"

class Acc{
	public:
	
	void toDist(std::string corrfname);
	double dist(double x1, double y1,double z1, double x2, double y2,double z2);	
	void dexus_cat_all(std::string outdir);
	void acc();
	void readDexus(std::string out_dir);
	void	readDexus_genes(std::string out_dir);
	std::vector<std::string> toPCA(std::string out_dir, int g_num);
		std::vector<std::string> toPCA_gene(std::string out_dir, int g_num);
	std::map<std::string, double> readAS(std::string allfname, std::string asfname);
	std::vector<std::string>	selectGene(std::vector<std::string> pathVec1,std::vector<std::string> pathVec2);
	std::vector<std::string>	readCuffdiff(std::string fname);
		
		
	void to_pROC_sdeap(std::string predfname, std::string allfname, std::string dsfname,std::string outfname);
	void to_PR_sdeap(std::string predfname, std::string allfname, std::string dsfname,std::string outfname);	
			
	void to_pROC_dexus(std::string predfname, std::string allfname, std::string dsfname,std::string outfname);	
	void to_PR_dexus(std::string predfname, std::string allfname, std::string dsfname,std::string outfname);	
	
	
		
};

#endif