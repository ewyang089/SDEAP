
#ifndef ASM_H
#define ASM_H
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
#include "SPGraph.h"

class ASM{
	public:
	std::map<std::string, SPGraph> ASM_Map;
	std::map<std::string, double> Junction_Map;	
	SPGraph Guide_Tr;
	void decomposeG(SPGraph InGraph);
	void evaluate_ASMs();	
	std::vector< std::vector<std::string> > toMTX_SPG(SPGraph spg, int len);
	
	private:
		
	std::vector<std::string>	_SetInOrder(std::vector<std::string> OrderVec, std::map<std::string, std::string> SetMap);
	
};

#endif