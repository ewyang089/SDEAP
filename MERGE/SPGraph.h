
#ifndef SPGRAPH_H
#define SPGRAPH_H
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
#include "Edge.h"
#include "Node.h"
#include "GenomeSegment.h"

class SPGraph{
	public:

	std::string gname; 
	int junction ;
	std::map<std::string, Node> NodeMap;
	std::map<std::string, Edge > EdgeMap;	   
		
		
  void _debug();
  void save(std::string FName);
  void load(std::string FName);
  	
  	
  int deleteEdge(std::string s_str, std::string t_str);
  int deleteNode(std::string s_str );			
  int addEdge(std::string s_str, std::string t_str, std::string edge_type,double weight ,int len );	
  int addNode(std::string s_str, double weight);	
  		
  std::vector<std::string> topological_sort();	
	void label_Pre_Dominator();
	void label_Post_Dominator();
	std::vector<std::string> list_Pred(std::string s_str);
	std::vector<std::string> list_Succ(std::string s_str);	
	
	std::vector<std::string> list_Paths();
	
	std::vector<std::string> list_UniPaths();
	
	
	SPGraph reverse(); 

	SPGraph subG(std::string s_str, std::string t_str); 
	SPGraph subASM(std::string s_str, std::string t_str); 
		
	SPGraph spanningT(std::string s_str); 
	
	
	std::vector<std::string> postorder(std::string s_str);
	std::vector< std::vector<std::string> > toMTX();
		
  void output_NodeWeight(std::string outfname);
	
	
	private:
	void _print_edge(std::string edgeneme);	
	std::map< std::string , std::string> map_intersect(std::map< std::string , std::string> AMap, std::map< std::string , std::string> BMap);		
	
};

#endif