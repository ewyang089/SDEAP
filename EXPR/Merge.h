#ifndef MERGE_H
#define MERGE_H
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

class Merge{
	public:
		
	std::vector< std::vector<double> > FeatureVec;
	std::vector< std::string > FNameVec;	
	std::vector< std::vector<double> > RCVec;
	std::vector< std::vector<double> > percentvec  ; 
		
		
	void merge4DexusExon(std::vector<std::string> pathvec, std::string outfname , std::string gene_list);
		
	void merge4DexusGene(std::vector<std::string> pathvec, std::string outfname , std::string gene_list);
	void merge_SPG(std::vector<std::string> pathvec, std::string outfname , std::string gene_list);	
	void merge_one_Hybrid(std::vector<std::string> pathvec, std::string outfname , std::string gene_list);	
	
	std::vector<std::string> mergeOneGene(std::vector<std::string> pathvec , std::string gene_list); 
	std::vector<std::string> mergeOneGene_select(std::vector<std::string> pathvec , std::string gene_list,double cut);
	std::vector<std::string> merge_Paths(std::vector<std::string> pathvec, std::string outdir , std::string gene_name,double cut);
	std::vector<std::string> merge_RCs(std::vector<std::string> pathvec, std::string outdir , std::string gene_name,double cut);	
	std::vector<std::string> mergeOnePath_select(std::vector<std::string> pathvec , std::string gene_name,double cut); 	 
	std::vector<std::string> mergeOneRC_select(std::vector<std::string> pathvec , std::string gene_name,double cut); 	 	
	std::vector<std::string> mergeOneGene4Dexus(std::vector<std::string> pathvec , std::string gene_list);	
	

	private:
			
};

#endif