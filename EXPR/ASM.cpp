#include "ASM.h"
using namespace std;

void ASM::decomposeG(SPGraph InGraph){
	SPGraph TmpGraph = InGraph;
	int ASM_count =0;
	vector<string> TVec = InGraph.topological_sort();
	TOK tok;
	stringstream ss;
	//find ASM
		
	
	for(int i =TVec.size()-1; i>=0 ; i-- ){
		if(TmpGraph.NodeMap.find(TVec[i])!= TmpGraph.NodeMap.end()){
			Node node_i = TmpGraph.NodeMap.find(TVec[i])->second;
		
			if(node_i.out_d > 1){
				for(int j = i+1 ; j<TVec.size(); j++){
					if(TmpGraph.NodeMap.find(TVec[j]) != TmpGraph.NodeMap.end()){
						Node node_j = TmpGraph.NodeMap.find(TVec[j])->second;
		
						if(node_j.in_d > 1){
					
					SPGraph ASM_G = TmpGraph.subG(node_i.name, node_j.name);
					
					if(ASM_G.NodeMap.size()>0){
					//check 3 asm properties 
						//cout  << "select" << node_i.name<< ","<< node_j.name<< endl;
					
						Node node_s  = ASM_G.NodeMap.find(node_i.name)->second;
						Node node_t  = ASM_G.NodeMap.find(node_j.name)->second;
					
					
						if(node_s.out_d > 1 && node_t.in_d >1){
						int all_in =1;
						map<string,Node>::iterator node_itr;
						for(node_itr = ASM_G.NodeMap.begin();node_itr != ASM_G.NodeMap.end();node_itr ++){
							Node node_u= TmpGraph.NodeMap.find(node_itr->first)->second;
							string node_name = node_itr->first;
							if(node_name.compare(node_s.name)!= 0 && node_name.compare(node_t.name)!= 0){
								map<string,string>::iterator str_itr;
								//check in_edges 
								for(str_itr= node_u.in_edgeMap.begin() ;str_itr != node_u.in_edgeMap.end() ; str_itr ++){
									vector<string> tokvec = tok.Tokenize(str_itr->first , "_");
									if(ASM_G.NodeMap.find(tokvec[0]) == ASM_G.NodeMap.end()){
										all_in = -1;
										//cout << str_itr->first << " invalid " << endl;
									}
								}
								// check out_edges 
								for(str_itr= node_u.out_edgeMap.begin() ;str_itr != node_u.out_edgeMap.end() ; str_itr ++){
									vector<string> tokvec = tok.Tokenize(str_itr->first , "_");
									if(ASM_G.NodeMap.find(tokvec[1]) == ASM_G.NodeMap.end()){
										all_in = -1;
										//cout << str_itr->first << " invalid " << endl;
									}
									
								}
								
								
							}
							
						}
						
						if(all_in ==1){
							//cout << node_s.name<<"," << node_t.name<< " is a ASM " << endl;
							
							ss.str("");
							ss<< node_s.name<< "_"<< node_t.name;
							if(ASM_G.EdgeMap.find(ss.str()) != ASM_G.EdgeMap.end()){
								ASM_G.junction =1;
								Edge &edge = ASM_G.EdgeMap.find(ss.str())->second;
								Junction_Map.insert(make_pair(node_s.name+ "~"+ node_t.name,edge.weight));
							}
							else{
								ASM_G.junction =0;
							}
							
							//create new node in the guided tree
							ss.str("");
							ss << node_s.name<<"~" << node_t.name;
							
							Guide_Tr.addNode(ss.str(),0);
							
							ASM_Map.insert(make_pair( ss.str(), ASM_G ));
							
							//chk containment
								
							map<string,Edge>::iterator edge_itr;
							for(edge_itr = ASM_G.EdgeMap.begin(); edge_itr != ASM_G.EdgeMap.end() ;edge_itr++ ){
								Edge asm_edge = edge_itr->second;
								if(asm_edge.type.compare("asm") == 0){
									//cout << node_s.name<<"," << node_t.name<< " contains " << asm_edge.name << endl;
									vector<string> tokvec = tok.Tokenize(asm_edge.name,"_");
									Guide_Tr.addEdge(ss.str(),tokvec[0]+"~"+tokvec[1],"",0,0);
								}
							}
							//shrink the s-t ASM into an edge;
							
													
							for(edge_itr = ASM_G.EdgeMap.begin() ;edge_itr != ASM_G.EdgeMap.end() ; edge_itr++){
								vector<string> tokvec = tok.Tokenize(edge_itr->first,"_");
								TmpGraph.deleteEdge(tokvec[0],tokvec[1]);
							}
														 
							for(node_itr = ASM_G.NodeMap.begin();node_itr != ASM_G.NodeMap.end();node_itr ++){
								string node_str = node_itr-> first;
								if(node_str.compare(node_s.name)!=0 && node_str.compare(node_t.name)!=0){
									TmpGraph.deleteNode(node_str);	
								}
							}	
							
							//calculate ASM edge weight
							
							double asm_weight =0;
							for(edge_itr = ASM_G.EdgeMap.begin() ;edge_itr != ASM_G.EdgeMap.end() ; edge_itr++){
								vector<string> tokvec = tok.Tokenize(edge_itr->first,"_");
								
								if(tokvec[0].compare(node_s.name)==0 && tokvec[1].compare(node_t.name)!=0){
									asm_weight +=(edge_itr->second).weight ;
								}
								if(tokvec[0].compare(node_s.name)!=0 && tokvec[1].compare(node_t.name)==0){
									asm_weight +=(edge_itr->second).weight ;
								}
								
							}
							asm_weight = asm_weight/2.0;
							TmpGraph.addEdge(node_s.name, node_t.name, "asm",asm_weight ,0 );	
							
							//
						}
						}
					}
						}
					}
				}
			}
		}
	}
	
	//cout << "ok ...." << endl;
	
	// finalize The guided tree by creating its root
	
	ss.str("");
	ss<<TVec[0] << "~"<< TVec[TVec.size()-1];
	Guide_Tr.addNode(ss.str(),0);
	map<string,Edge>::iterator edge_itr;
  for(edge_itr = TmpGraph.EdgeMap.begin(); edge_itr != TmpGraph.EdgeMap.end() ;edge_itr++ ){
		Edge asm_edge = edge_itr->second;
		if(asm_edge.type.compare("asm") == 0){
			
			vector<string> tokvec = tok.Tokenize(asm_edge.name,"_");
			Guide_Tr.addEdge(ss.str(),tokvec[0]+"~"+tokvec[1],"",0,0);
		}
	}
	
	ASM_Map.insert(make_pair( ss.str(), TmpGraph ));
	
}
void ASM::evaluate_ASMs(){
	vector<string> postvec = Guide_Tr.postorder("S~T");
	for(int i=0 ; i< postvec.size();i++){
		map<string, SPGraph>::iterator asm_itr = ASM_Map.find(postvec[i]); 
		if( asm_itr != ASM_Map.end()){
			SPGraph spg = asm_itr-> second;
			vector<vector<string> > mtxvec = spg.toMTX();
			for(int j = 0 ; j < mtxvec.size() ; j++){
				for(int k =0; k< mtxvec[j].size(); k++){
					cout << mtxvec[j][k] << ",";
				}
				cout << endl;
			}
		}
	}
}

std::vector<std::vector<string> > ASM::toMTX_SPG(SPGraph spg, int len){
	//search the start && end node;

	map<string, Node>::iterator node_itr;
	string s_str,t_str;
	TOK tok;
	std::vector< vector<string> > outVec;	
	
	stringstream ss;
	
	
	for(node_itr = spg.NodeMap.begin(); node_itr != spg.NodeMap.end(); node_itr++){
		
		Node & node_u = node_itr->second;
		if(node_u.in_d ==0){
			s_str= node_itr->first;
		}
		if(node_u.out_d == 0){
			t_str= node_itr->first;
		}
	}
	



	//cout << "S:"<< s_str<< endl;
	//cout << "T:"<< t_str<< endl;	 
	vector<string> pathvec = spg.list_Paths();
/*
for(int i =0;i<pathvec.size() ; i++){
	cout <<"path:" <<pathvec[i]<< endl; 
}
*/

	list<string> Q;
	
	for(int i =0 ; i<pathvec.size();i++){
		vector<string> tokvec = tok.Tokenize(pathvec[i],",");
		ss.str("");
		if(tokvec.size()>=2){
			for(int j =0 ; j<tokvec.size(); j++){
				
					Node &node_i = spg.NodeMap.find(tokvec[j])->second ; 
			
					ss<< "N:"<< node_i.name << ":"<<node_i.weight;
			
					if(j<tokvec.size()-1){
				
					Edge &edge_e = spg.EdgeMap.find(tokvec[j]+"_"+tokvec[j+1])->second ; 
				
					if(edge_e.type.compare("asm") == 0){
					
					//ss<< ",asm:"<< tokvec[j]<< "~"<< tokvec[j+1] << ":"<< edge_e.weight<<",";
					ss<<",";
					}
					else{
					// edge_e
					ss<< ",J:"<< edge_e.name << ":"<<edge_e.weight<<",";
					
					}
				
						
					}
			
			}
				Q.push_back(ss.str());
		}
	}
	vector<string> SavePathVec;
	
	for(list<string>::iterator Q_itr = Q.begin() ;Q_itr != Q.end(); Q_itr++ ){
		//cout << "Q:"<<*Q_itr<<endl;
		SavePathVec.push_back(*Q_itr);
	}

	
	//parse segs
	
	
	map<string,string> segmap;
	map<string,string> Rmap;
	
		//set col names
	vector<string> col_name_vec;
	col_name_vec.push_back(" ");
	
	for(int i =0 ; i<SavePathVec.size() ;i++){
		vector<string> tokvec = tok.Tokenize(SavePathVec[i],",");
		
		col_name_vec.push_back(SavePathVec[i]);
		
		for(int j =1 ; j< tokvec.size()-1 ; j++){
			
			if(segmap.find(tokvec[j])== segmap.end()){
				segmap.insert(make_pair( tokvec[j],""));
			}
			vector<string> rvec  = tok.Tokenize(tokvec[j],":");
			if(Rmap.find(tokvec[j])== Rmap.end()){
				Rmap.insert(make_pair(tokvec[j],rvec[2]));
			}
		}
						
	}
	col_name_vec.push_back("R");
	outVec.push_back(col_name_vec);
	
	//create init rows
	map<string, int> SegIdxMap;
	int idx=1;
	for(map<string,string>::iterator seg_itr = segmap.begin(); seg_itr != segmap.end(); seg_itr++){
		vector<string> tmpvec ;
		tmpvec.push_back(seg_itr->first);
		for(int i =0 ; i<SavePathVec.size() ;i++){
			tmpvec.push_back("0");
		}
		tmpvec.push_back(Rmap.find(seg_itr->first)->second);
		outVec.push_back(tmpvec);
		SegIdxMap.insert( make_pair(seg_itr->first , idx) );
		idx++;
		}
	
	//set length in mtx
		for(int i =0 ; i<SavePathVec.size() ;i++){
			vector<string> tokvec = tok.Tokenize(SavePathVec[i],",");
		
			for(int j =1 ; j< tokvec.size()-1 ; j++){
			
			vector<string> rvec = tok.Tokenize(tokvec[j],":");
			int idx = SegIdxMap.find(tokvec[j])->second;
				ss.str("");
				if(rvec[0].compare("N")==0){
					vector<string> posvec = tok.Tokenize(rvec[1],"$");
				//cout << "pos 1:"<< posvec[0]<< endl;
			 //cout << "pos 2:"<< posvec[1]<< endl;
					ss<<(atoi(posvec[1].c_str())-atoi(posvec[0].c_str())+1)/(double)1000.0 ;
				}
				else{
					ss<<len/(double)1000;
				}
				outVec[idx][i+1] = ss.str();
			
			}
						
		}
	
	
	
	return outVec;
}

std::vector<string>	ASM::_SetInOrder(std::vector<string> OrderVec, std::map<std::string, std::string> SetMap){
	
}