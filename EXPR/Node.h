
#ifndef NODE_H
#define NODE_H
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

class Node{
	public:

	std::string name;
	std::string type;
	double weight;	
	int in_d;
	int out_d;  	
	std::map<std::string, std::string> out_edgeMap,in_edgeMap;
	std::map<std::string, std::string> PreDomMap,PostDomMap;
			
	private:
			
};

#endif