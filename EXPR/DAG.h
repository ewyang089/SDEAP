
#ifndef DAG_H
#define DAG_H
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
#include "factor.h"

class DAG{
	public:
  DAG(); 
  std::map<std::string, double > nodeMap;
	std::map<std::string, std::map<std::string, double> > edgeMap;	
	std::map<std::string, factor > bpMap;
	DAG merge(DAG ADAG,DAG BDAG);	
	void write(std::string FName);
	void read(std::string FName);
	std::vector<std::string> mincut(std::string sstr, std::string tstr);	
	void print();
	double getPOSPROB(std::string nname);
	std::vector<DAG> split();
	void bp(std::string SourceNode,int round);
	void bp_revised(std::string SourceNode,int round);
	void bp_noiseq(std::string SourceNode,int round);
	std::vector<std::string> scheduleNodes(std::string SourceNode);
	
	private:
	std::map<std::string, int> GtoXMap;
	std::vector<std::string> partition(std::string source,DAG &newdag);
	DAG antiparaedgefreeDAG();
	DAG residueG();	
	std::vector<std::string> BFS(std::string sstr, std::string tstr, DAG &inputG);
	void convertFactors();
		
};

#endif
