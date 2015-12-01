#include "SPGraph.h"
using namespace std;

void SPGraph::label_Post_Dominator(){
	SPGraph spg = reverse();
	spg.label_Pre_Dominator();
	map<string, Node>::iterator node_itr ;
	for(node_itr = spg.NodeMap.begin() ;node_itr != spg.NodeMap.end() ;node_itr ++){
		Node & node_t = NodeMap.find(node_itr->first)->second;
		node_t.PostDomMap = (node_itr->second).PreDomMap;
	}
	
}


	
	
std::vector<string> SPGraph::list_Paths(){
	vector<string> OutVec;
	TOK tok;
  map<string, vector<string> > PathMap;
  vector<string> NodeNameVec= topological_sort();
   
  //init  PathMap 
  for(int i =0; i < NodeNameVec.size(); i++){
  	vector<string> tmpvec ;
  	PathMap.insert(make_pair(NodeNameVec[i],tmpvec));
  }
   	
  map<string, vector<string> >::iterator node_itr = PathMap.find(NodeNameVec[0]);
  (node_itr->second).push_back(NodeNameVec[0]);
  
  for(int i =0; i < NodeNameVec.size(); i++){
  	map<string, vector<string> >::iterator u_itr = PathMap.find(NodeNameVec[i]);
  		Node & node_u = NodeMap.find(NodeNameVec[i])->second;
  		for(map<string,string>::iterator edge_itr =node_u.in_edgeMap.begin()  ;edge_itr !=node_u.in_edgeMap.end() ; edge_itr++){
			
				vector<string> tokvec = tok.Tokenize(edge_itr->first,"_");
				vector<string> pvec = PathMap.find(tokvec[0])->second;
				for(int j=0;j<pvec.size();j++){
					(u_itr->second).push_back(pvec[j]+","+NodeNameVec[i]);
				}
			}	
  		
  }
  OutVec = PathMap.find(NodeNameVec[NodeNameVec.size()-1])->second	; 		
	return OutVec;
	
}

std::vector<string> SPGraph::list_UniPaths(){

	TOK tok;
  stringstream ss;
  vector<string> PathVec;
  
  map<string, string> startmap;
  
  //label nodes with degree >1
  for(map<string, Node>::iterator node_itr = NodeMap.begin() ; node_itr != NodeMap.end(); node_itr++){
  	Node &node_u = node_itr->second;
  	if( node_u.in_d ==0 || node_u.in_d > 1  || node_u.out_d > 1){
  		startmap.insert(make_pair(node_itr->first,""));
  	}
  	
  }
  
  
  
	for( map<string, string>::iterator str_itr = startmap.begin()  ; str_itr != startmap.end(); str_itr++){
		SPGraph T = spanningT(str_itr->first);
		//cout << "####################################"<< endl;
		
		//DFS to find all the paths from the root		
		vector<string> S ; 
		map<string,string> addedmap;
		S.push_back(str_itr->first);
		addedmap.insert(make_pair(str_itr->first ,""));
		
		while(S.size() > 0 ){
			string u_str = S[S.size()-1];
			
			Node &node_u = T.NodeMap.find(u_str)->second;
			Node &_u = NodeMap.find(u_str)->second;
			int in_d = _u.in_d;
			int out_d =_u.out_d ;
			//cout <<u_str <<"," << in_d<< ","<< out_d << "," << str_itr->first << endl;
			if( ((in_d > 1|| out_d >1) && u_str.compare(str_itr->first) != 0) || node_u.out_d == 0  ){
				
				ss.str("");
				for(int i  = 0; i < S.size()  ; i++){
						ss<< S[i];
						if(i<S.size()-1){
							ss<<",";
						}
					//cout << "add path:" <<S[i] << endl;
				}
				S.pop_back();
				PathVec.push_back(ss.str());
				
			}
			else{
				int flag =0;
				for(map<string,string>::iterator v_itr = node_u.out_edgeMap.begin();v_itr != node_u.out_edgeMap.end() ;v_itr++ ){
					vector<string> posvec= tok.Tokenize(v_itr->first,"_");
					
					if(addedmap.find(posvec[1]) == addedmap.end()){
						S.push_back(posvec[1]);
						addedmap.insert(make_pair(posvec[1],""));
						
						flag =1;
						break;
					}
				}
				
				if(flag ==0 ){
					S.pop_back();
				}
			}
			
			
		}
		
	}
	
	return PathVec;
	
}

SPGraph SPGraph::subASM(std::string s_str, std::string t_str){
	SPGraph outGraph; 
	
	
	if(s_str.compare(t_str) ==0){
		return outGraph;
	}
	if(NodeMap.find(s_str) == NodeMap.end()){
		return outGraph;
	}
	if(NodeMap.find(t_str) == NodeMap.end()){
		return outGraph;
	}
	Node node_s = NodeMap.find(s_str)->second;
	if(node_s.PostDomMap.find(t_str) == node_s.PostDomMap.end()){
		return outGraph;
	}
	
	Node node_t = NodeMap.find(t_str)->second;
	if(node_t.PreDomMap.find(s_str) == node_t.PreDomMap.end()){
		return outGraph;
	}
	
	
	TOK tok;
	vector<string> S;
	map<string , string > AddedMap ; 
	
	S.push_back(s_str);
	AddedMap.insert(make_pair( s_str, ""));
  AddedMap.insert(make_pair( t_str, ""));
	while(S.size()>0){
		string u_str = S[0];
		S.erase(S.begin());
		Node node_u = NodeMap.find(u_str)->second;
		if(u_str.compare(s_str)==0){
			node_u.in_edgeMap.clear();
			node_u.in_d = 0;
		}
		
		outGraph.NodeMap.insert(make_pair( u_str,node_u ));
		
		map<string,string>::iterator out_itr;
		
		for(out_itr = node_u.out_edgeMap.begin();out_itr !=node_u.out_edgeMap.end() ; out_itr ++){
			string outedge_str = out_itr -> first;
			
			Edge outedge  = (EdgeMap.find(outedge_str))->second;

			outGraph.EdgeMap.insert(make_pair( outedge_str , outedge ));

			vector<string> tokVec = tok.Tokenize(outedge_str,"_");

			if(AddedMap.find(tokVec[1]) ==AddedMap.end() ){

				AddedMap.insert(make_pair( tokVec[1], ""));
				S.push_back(tokVec[1] );

			}
		}
		
	}
	Node Node_t = NodeMap.find(t_str)->second;
	Node_t.out_edgeMap.clear();
	Node_t.out_d = 0;
	outGraph.NodeMap.insert(make_pair( t_str, Node_t));
	return outGraph;
	
} 


SPGraph SPGraph::spanningT(std::string s_str){
	SPGraph G;	
	vector<string> List;
	TOK tok ;
	map<string, string> AddedMap; 
	stringstream ss;
	List.push_back(s_str);
		while(List.size()>0){
	
	
		string u_str = List[0];
		List.erase(List.begin());
	  
	  AddedMap.insert(make_pair(s_str,""));
	  Node & node_u = NodeMap.find(u_str)->second;
		G.addNode(u_str, node_u.weight);
		
		map<string , string>::iterator outedge_itr;

		for(outedge_itr = node_u.out_edgeMap.begin(); outedge_itr != node_u.out_edgeMap.end(); outedge_itr++){
			
			Edge edge_uv = (EdgeMap.find(outedge_itr->first))->second;
			vector<string> tokvec = tok.Tokenize( outedge_itr->first, "_");
			
			if(AddedMap.find(tokvec[1]) == AddedMap.end()){
				// add the edge to the tree
				Node & node_v = NodeMap.find(tokvec[1])->second; 
				G.addNode(tokvec[1], node_v.weight);
				G.addEdge(tokvec[0], tokvec[1], edge_uv.type, edge_uv.weight ,edge_uv.len );	
				List.push_back(tokvec[1]);
				AddedMap.insert(make_pair(tokvec[1],""));
				
			}
			
		} 
		
	}
	return G;
}

std::vector<std::vector<string> > SPGraph::toMTX(){
	
	//search the start && end node;
	map<string, Node>::iterator node_itr;
	string s_str,t_str;
	TOK tok;
	std::vector< vector<string> > outVec;	
	stringstream ss;
	
	
	for(node_itr = NodeMap.begin(); node_itr != NodeMap.end(); node_itr++){
		
		Node & node_u = node_itr->second;
		if(node_u.in_d ==0){
			s_str= node_itr->first;
		}
		if(node_u.out_d == 0){
			t_str= node_itr->first;
		}
	}
	


	//find how many segs
	vector<string> pathvec = list_Paths();
	map<string,string> segmap;

	
	int col_n = 0;
	int row_n = 0 ;
	
	//set col names
	vector<string> col_name_vec;
	col_name_vec.push_back(" ");
	
	for(int i =0 ; i<pathvec.size() ;i++){
		
		vector<string> tmpvec = tok.Tokenize(pathvec[i],",");
		if(tmpvec.size()>2){
			col_n++;
			col_name_vec.push_back(pathvec[i]);
			for(int j =1 ; j< tmpvec.size()-1 ; j++){
				if(segmap.find(tmpvec[j])== segmap.end()){
					segmap.insert(make_pair( tmpvec[j],""));
				}
			}
						
		}
		
	}
	
	col_name_vec.push_back("R");
  outVec.push_back(col_name_vec);
	
	row_n = segmap.size();
	row_n = row_n +1;

  
	map<string, int > seg_idx_map;
	int srg_idx =1;  
  
  for(map<string,string>::iterator seg_itr = segmap.begin() ;seg_itr != segmap.end() ;seg_itr++ ){
  	vector<string> segvec ;
  	segvec.push_back(seg_itr->first);
  	seg_idx_map.insert(make_pair(seg_itr->first, srg_idx ));
  	for(int i=0;i<col_n;i++){
			segvec.push_back("0");
			
		}
  	segvec.push_back("0");
  	outVec.push_back(segvec);
  	srg_idx++;
  }
  
  //fill length of the seg
  int path_idx =1;
  for(int j =0 ; j < pathvec.size() ; j++){
  	vector<string> tmpvec = tok.Tokenize(pathvec[j],",");
		if(tmpvec.size()>2){
			
			for(int q =0 ; q< tmpvec.size() ; q++){
				map<string, int >::iterator segidx_itr =  seg_idx_map.find(tmpvec[q]);
				if(segidx_itr!= seg_idx_map.end()){
					int i = segidx_itr->second;
					vector<string> posvec = tok.Tokenize(tmpvec[q] , "$");
					int spos = atoi(posvec[0].c_str());
					int epos = atoi(posvec[1].c_str());
					ss.str("");
					
					ss<<epos-spos+1;
					outVec[i][path_idx] = ss.str();
				}	
			}
			path_idx++;			
		}
  	
  }
  
  //fill observed R
  
  for(map<string, int >::iterator segidx_itr =  seg_idx_map.begin(); segidx_itr !=  seg_idx_map.end(); segidx_itr++){
  	
  	map<string, Node>::iterator node_itr = NodeMap.find(segidx_itr->first);
  	if(node_itr != NodeMap.end()){
  		Node & node_u  = node_itr->second;
  	  ss.str("");
  	  ss<< node_u.weight;
  	  outVec[segidx_itr->second][path_idx] = ss.str();
  	}
  
  }
  
  
	return outVec;
}

vector<string> SPGraph::postorder(std::string s_str){
	vector<string> outVec; 
	TOK tok;
	list<string> S;
	S.push_front(s_str);
	map<string,string> AddedMap;
		
	while(S.size()>0){
		string u_str = S.front();
		
		Node node_u = NodeMap.find(u_str)->second;
		AddedMap.insert(make_pair(u_str,""));
		
		int flag =1 ;
		
		map<string,string>::iterator node_itr ;
		
		for(node_itr = node_u.out_edgeMap.begin();node_itr != node_u.out_edgeMap.end();node_itr++){
			vector<string> tokvec = tok.Tokenize(node_itr->first ,"_");
			if(AddedMap.find(tokvec[1]) == AddedMap.end()){
				
				S.push_front(tokvec[1]);
				flag =-1;
				break;
			}
			
		}
		
		if(flag ==1){
			outVec.push_back(u_str);
			
			S.pop_front();
		}
		
	}
	
	return outVec;
}
	
	
SPGraph SPGraph::subG(std::string s_str, std::string t_str){
	
	
	//the Spaning tree from s_str
	TOK tok ;
	SPGraph sT = spanningT(s_str);
	
	SPGraph Reversed_G = reverse(); 
	
	SPGraph tT =Reversed_G.spanningT(t_str);
	
	if(sT.NodeMap.find(t_str)==sT.NodeMap.end()){
		SPGraph emptyG;
		return emptyG;
	}
	
	
	map<string , Node>::iterator node_itr; 
	SPGraph OutGraph ;
	
	for(node_itr = sT.NodeMap.begin();node_itr != sT.NodeMap.end(); node_itr ++){
		if(tT.NodeMap.find(node_itr->first)!= tT.NodeMap.end()){
			Node & node_u = node_itr -> second; 
			OutGraph.addNode(node_u.name,node_u.weight);
		}
			
	}
	
	 
	//pick edges
	
	vector<string> nodevec ;
	for(node_itr = OutGraph.NodeMap.begin(); node_itr != OutGraph.NodeMap.end(); node_itr++){
		nodevec.push_back(node_itr->first);
	}
			
	for(int i =0 ; i<nodevec.size();i++ ){
		Node node_u = NodeMap.find(nodevec[i])->second;
		//add in_edges
		
		map<string, string>::iterator edge_itr;
		for(edge_itr = node_u.in_edgeMap.begin() ; edge_itr != node_u.in_edgeMap.end() ; edge_itr++){
		
			vector<string> tokvec = tok.Tokenize( edge_itr->first, "_" );
			if( OutGraph.NodeMap.find(tokvec[0]) !=OutGraph.NodeMap.end() && OutGraph.NodeMap.find(tokvec[1]) !=OutGraph.NodeMap.end()){
				
				Edge & newedge = EdgeMap.find(edge_itr->first)->second;
				
				if(OutGraph.EdgeMap.find(edge_itr->first) == OutGraph.EdgeMap.end()){
				OutGraph.addEdge(tokvec[0], tokvec[1], newedge.type,newedge.weight ,newedge.len );
				
				}
			}
			
		}
			
		
		//add out_edges
		for(edge_itr = node_u.out_edgeMap.begin() ; edge_itr != node_u.out_edgeMap.end() ; edge_itr++){
			vector<string> tokvec = tok.Tokenize( edge_itr->first, "_" );
			if( OutGraph.NodeMap.find(tokvec[0]) !=OutGraph.NodeMap.end() && OutGraph.NodeMap.find(tokvec[1]) !=OutGraph.NodeMap.end()){
				Edge & newedge = EdgeMap.find(edge_itr->first)->second;
				if(OutGraph.EdgeMap.find(edge_itr->first) == OutGraph.EdgeMap.end()){			
					OutGraph.addEdge(tokvec[0], tokvec[1], newedge.type,newedge.weight ,newedge.len );
				}
			}
		}
			
	}
	
	return OutGraph;
	
	
}

SPGraph SPGraph::reverse(){
	SPGraph OutGraph;
	TOK tok;
	stringstream ss;
	map<string , Node>::iterator nodemap_itr ;
		 
	for(nodemap_itr =  NodeMap.begin();nodemap_itr !=  NodeMap.end() ; nodemap_itr++){
		Node newNode ;
		Node &Node_t = nodemap_itr->second;
		 
		newNode.name =Node_t.name ;
		newNode.type =Node_t.type ;
		
		newNode.weight =Node_t.weight ;
		newNode.in_d =Node_t.out_d ;
		newNode.out_d =Node_t.in_d ;
		
		map<string, string>::iterator mapitr ; 
	
		for(mapitr = Node_t.in_edgeMap.begin() ;  mapitr != Node_t.in_edgeMap.end(); mapitr ++){
			vector<string> tokvec = tok.Tokenize(mapitr->first,"_");
			ss.str("");
			ss << tokvec[1]<< "_"<< tokvec[0];
			newNode.out_edgeMap.insert(make_pair( ss.str(), ""));
		}
	
		
		for(mapitr = Node_t.out_edgeMap.begin() ;  mapitr != Node_t.out_edgeMap.end(); mapitr ++){
			vector<string> tokvec = tok.Tokenize(mapitr->first,"_");
			ss.str("");
			ss << tokvec[1]<< "_"<< tokvec[0];
			newNode.in_edgeMap.insert(make_pair( ss.str(), ""));
		}

		OutGraph.NodeMap.insert(make_pair( newNode.name , newNode));	
	}
	
	
	map<string , Edge>::iterator edgemap_itr ;
		
	for(edgemap_itr = EdgeMap.begin() ; edgemap_itr != EdgeMap.end(); edgemap_itr ++){
		Edge newEdge ;
		Edge &Edge_t = edgemap_itr->second;
		
		vector<string> tokvec = tok.Tokenize(Edge_t.name,"_");
		ss.str("");
		ss << tokvec[1]<< "_"<< tokvec[0];
		newEdge.name = ss.str();
		newEdge.type = Edge_t.type;
		newEdge.len = Edge_t.len;
		newEdge.weight = Edge_t.weight;
		
		OutGraph.EdgeMap.insert(make_pair( newEdge.name , newEdge));	
	}
	
	return OutGraph;
}
void SPGraph::_debug(){

	string str;	 
	vector<string> tokvec;
	TOK tok;
	map<string,Node>::iterator nodeitr;
	map<string, Edge >::iterator edgeitr;	

  cout << "SIZE:" << NodeMap.size()<< ","<<EdgeMap.size() <<  endl;
	cout << "Jucntoion:" << junction <<  endl;
	for(nodeitr = NodeMap.begin();nodeitr != NodeMap.end() ; nodeitr++){
		Node & node = nodeitr->second;
		cout << endl<<"NODE :  " << node.name <<  ":"<<node.weight<<":"<<node.in_d << ":" << node.out_d << endl;
		
		cout << "PreDom :  "<<endl;
		
		map<string,string>::iterator predom_itr;
		
		for( predom_itr = node.PreDomMap.begin();predom_itr != node.PreDomMap.end() ; predom_itr++){
			cout << predom_itr -> first <<endl;
			
		}
		
		cout << "PostDom :  "<<endl;
		
		map<string,string>::iterator postdom_itr;
		
		for( postdom_itr = node.PostDomMap.begin();postdom_itr != node.PostDomMap.end() ; postdom_itr++){
			cout << postdom_itr -> first <<endl;
			
		}
		
		map<string,string>::iterator edge_itr;
		
		
		cout << "In_Edges :  "<<endl;
		
		for(edge_itr =node.in_edgeMap.begin()  ;edge_itr !=node.in_edgeMap.end() ; edge_itr++){
			
			_print_edge(edge_itr->first);
		}	
		
		cout << "Out_Edges :  "<<endl;
		
		for(edge_itr =node.out_edgeMap.begin()  ;edge_itr !=node.out_edgeMap.end() ; edge_itr++){
			_print_edge(edge_itr->first);
		}
	}
	
}

void SPGraph::_print_edge(string edgeneme){
	Edge & edge = EdgeMap.find(edgeneme)->second;
	cout <<edge.name << ":"<< edge.type << ":"<< edge.weight << ":"<< edge.len<<endl;
}

int SPGraph::addNode(std::string s_str, double weight  ){
	int outInt =-1;
	
	Node NewNode ; 	
	
	if( NodeMap.find(s_str) == NodeMap.end() ){
		NewNode.name = s_str;
		NewNode.weight = weight;
		NewNode.in_d = 0;
		NewNode.out_d = 0;
		outInt =1;
		NodeMap.insert(make_pair(s_str, NewNode));
		
	}
	return outInt;
}

std::vector<std::string> SPGraph::list_Pred(std::string s_str){
	vector<string> outVec;
	TOK tok;
	Node & node = NodeMap.find(s_str)->second;
	map<string, string> in_edgeMap = node.in_edgeMap;
	map<string, string>::iterator edgeitr;
		for(edgeitr = in_edgeMap.begin(); edgeitr != in_edgeMap.end(); edgeitr++){
			vector<string> tokvec = tok.Tokenize(edgeitr->first,"_");
			
			outVec.push_back(tokvec[0]);
		}
		
		
	return outVec;
}


std::vector<std::string> SPGraph::list_Succ(std::string s_str){
	vector<string> outVec;
	TOK tok;
	Node & node = NodeMap.find(s_str)->second;
	map<string, string> out_edgeMap = node.out_edgeMap;
	map<string, string>::iterator edgeitr;
		for(edgeitr = out_edgeMap.begin(); edgeitr != out_edgeMap.end(); edgeitr++){
			vector<string> tokvec = tok.Tokenize(edgeitr->first,"_");
			
			outVec.push_back(tokvec[1]);
		}
		
		
	return outVec;
}

std::vector<string> SPGraph::topological_sort(){
	vector<string> S,L ;
	
	map<string,Node>::iterator nodeitr;
	//initialize S ;	
	for(nodeitr=NodeMap.begin() ; nodeitr!=NodeMap.end(); nodeitr++){
		Node & node = nodeitr -> second ;  
		if(node.in_d == 0){
			S.push_back(node.name);
			
		}
	}
	
	map<string, string> visitedMap;
	
	while(S.size()>0){
		string node_str = S[0];
		L.push_back(S[0]);
		S.erase(S.begin());
		visitedMap.insert(make_pair(node_str,""));		  
	
		vector<string> SuccVec = list_Succ(node_str);
		
		for(int i =0;i<SuccVec.size(); i++){
			string t_str = SuccVec[i];
			vector<string> Pred = list_Pred(t_str);
			int flag = -1;
			for(int j=0 ; j<Pred.size();j++){
				if(visitedMap.find(Pred[j]) == visitedMap.end()){
					flag =1;
				}
			}
			
			if(flag == -1){
				S.push_back(t_str);
			}
		} 
		
	}
	
	
	return L ;
}	

void SPGraph::label_Pre_Dominator(){
	vector<string> OrderVec = topological_sort();
	string  str_s;
	//init 
	for(int i =0; i < OrderVec.size(); i++){
		 str_s = OrderVec[i];
		//Node & node = NodeMap.find(node_s)->second;
		cout << "label :" << str_s<< endl;
		vector<string> PredVec = list_Pred(str_s);
		
		map<string, string> PreDomMap;
		
			for(int j = 0 ; j < PredVec.size() ; j++){
				
			Node & node_t = NodeMap.find(PredVec[j])->second;
			
			if(j==0){
				PreDomMap = (NodeMap.find(PredVec[j])->second).PreDomMap;
			}
			else{
			  PreDomMap = map_intersect(PreDomMap, (NodeMap.find(PredVec[j])->second).PreDomMap);
			}
		}		
		
		PreDomMap.insert( make_pair( str_s, "") );
	  //Node & node_s = NodeMap.find(str_s)->second;
	  (NodeMap.find(str_s)->second).PreDomMap = PreDomMap;
	}
	
}

std::map< std::string , std::string> SPGraph::map_intersect(std::map< std::string , std::string> AMap, std::map< std::string , std::string> BMap){
	map<string, string> OutMap;
	map<string, string>::iterator map_itr; 
	cout << "inter ?? " <<AMap.size()<< "" <<BMap.size() << endl;
	for(map_itr = AMap.begin(); map_itr != AMap.end(); map_itr++){
		if(BMap.find(map_itr->first)!= BMap.end()){
			OutMap.insert(make_pair( map_itr->first , ""));
		}
	}
	
	return OutMap; 
	
}		

int SPGraph::deleteNode(std::string s_str){
	int outInt =-1;
	TOK tok;
	if(NodeMap.find( s_str ) != NodeMap.end()){
		Node node = NodeMap.find( s_str )->second;
		map<string,string>::iterator edgeitr ;
			
		for(edgeitr = node.out_edgeMap.begin();edgeitr != node.out_edgeMap.end(); edgeitr ++){
			string edge_name = edgeitr->first;
			vector<string> tokvec = tok.Tokenize(edge_name, "_");
			deleteEdge(tokvec[0], tokvec[1] );
		}
		
		for(edgeitr = node.in_edgeMap.begin();edgeitr != node.in_edgeMap.end(); edgeitr ++){
			string edge_name = edgeitr->first;
			vector<string> tokvec = tok.Tokenize(edge_name, "_");
			deleteEdge(tokvec[0], tokvec[1] );
		}
		
		outInt =1;
		NodeMap.erase(s_str);
	}
	
	return outInt;
}

int SPGraph::deleteEdge(std::string s_str, std::string t_str){
	int outInt =-1;
	
	stringstream ss;
	ss << s_str<< "_"<< t_str ;
	string edgename = ss.str(); 
	if(EdgeMap.find( edgename ) != EdgeMap.end()){
		Node  node = NodeMap.find(s_str)->second;
		node.out_d--;
		node.out_edgeMap.erase(edgename);
		NodeMap.erase(s_str);
		NodeMap.insert( make_pair( s_str, node) );
		
		
		node = NodeMap.find(t_str)->second;
		node.in_d--;
		node.in_edgeMap.erase(edgename);
		NodeMap.erase(t_str);
		NodeMap.insert( make_pair( t_str, node) );
		EdgeMap.erase(edgename);
	
		outInt =1;	
	}
	
	return outInt;
}

int SPGraph::addEdge(std::string s_str, std::string t_str,  std::string edge_type,double weight ,int len ){
	int outInt =-1;
	Edge NewEdge ; 
	stringstream ss;
	
	
	if( NodeMap.find(s_str) != NodeMap.end() && NodeMap.find(t_str) != NodeMap.end()){
		ss.str("");
		ss<< s_str << "_"<< t_str;
		if(EdgeMap.find(ss.str())==EdgeMap.end() ){
			outInt =1;
			NewEdge.name = ss.str();
			NewEdge.type = edge_type;
			NewEdge.weight = weight;
			NewEdge.len = len;
			EdgeMap.insert(make_pair(ss.str(), NewEdge));
		
			Node node_s = NodeMap.find(s_str)->second;
			node_s.out_d++;
			node_s.out_edgeMap.insert(make_pair(ss.str(), ""));
		
			NodeMap.erase(s_str);
			NodeMap.insert(make_pair(s_str,node_s));
		
			Node node_t = NodeMap.find(t_str)->second;
			node_t.in_d++;		
			node_t.in_edgeMap.insert(make_pair(ss.str(), ""));
			NodeMap.erase(t_str);
			NodeMap.insert(make_pair(t_str,node_t));
		}
	}
	return outInt;
}

void SPGraph::save(string FName){
	fstream filestr;
	filestr.open(FName.c_str(), ios_base::out);
	string str;	 
	vector<string> tokvec;

	map<string, Node>::iterator node_itr;
	map<string, Edge>::iterator edge_itr;
	for(node_itr = NodeMap.begin();node_itr != NodeMap.end();node_itr++){
		Node &node_u = node_itr->second;
		filestr <<  "N:"<< node_u.name<< ","<<node_u.weight<<endl;
	}	
	
	for(edge_itr = EdgeMap.begin() ;edge_itr != EdgeMap.end(); edge_itr++){
		Edge & edge_e = edge_itr->second;
		filestr <<  "E:"<< edge_e.name<< ","<< edge_e.type<< "," <<edge_e.weight<< ","<< edge_e.len <<endl;  
	}
	filestr.close();
}

void SPGraph::output_NodeWeight(std::string outfname )
{
	map<string, Node>::iterator node_itr;
		TextIO tio;
		stringstream ss;
		vector<string> outVec;
		for(node_itr = NodeMap.begin(); node_itr != NodeMap.end(); node_itr++){
			Node & node_u = node_itr->second;
			ss.str("");
			ss<<gname<<":"<< node_itr->first <<"\t" <<node_u.weight;
			outVec.push_back(ss.str() );
		}
		tio.dlwrite(outfname, outVec);
		 
}
void SPGraph::load(string FName){
	fstream filestr;
	string str;	 
	vector<vector<string> > tokvec;
	TextIO tio;
	TOK tok;
	tokvec = tio.dlread( FName,":");
	
	for(int i =0; i<tokvec.size(); i++){
		
		if(tokvec[i][0].compare("N")==0){
			vector<string> tmpvec = tok.Tokenize(tokvec[i][1],",");
			
			addNode(tmpvec[0], atof(tmpvec[1].c_str()));
		}
		if(tokvec[i][0].compare("E")==0){
			vector<string> tmpvec = tok.Tokenize(tokvec[i][1],",");
			vector<string> nodevec= tok.Tokenize(tmpvec[0],"_");
			addEdge(nodevec[0], nodevec[1],tmpvec[1],atof(tmpvec[2].c_str()), atoi(tmpvec[3].c_str()) );
		}
		
	}
	
	
}
