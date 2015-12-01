
#include "DAG.h"

using namespace std;

DAG::DAG(){
	
}

void DAG::convertFactors(){
	
}

void DAG::write(string FName){
	fstream filestr;
	filestr.open(FName.c_str(), ios_base::out);
	string str;	 
	vector<string> tokvec;
	TOK tok;
	map<string,double>::iterator nodeitr;
	map<string, map<string, double> >::iterator edgeitr;	
	for(nodeitr = nodeMap.begin();nodeitr != nodeMap.end() ; nodeitr++){
		filestr<< nodeitr->first<<"_"<<nodeitr->second << ":";
		edgeitr = edgeMap.find(nodeitr->first);
		map<string,double> LinkMap = edgeitr->second;
		map<string,double>::iterator LinkItr;
		for(LinkItr = LinkMap.begin(); LinkItr != LinkMap.end(); LinkItr++){
			filestr<< LinkItr->first << "_"<< LinkItr->second << ",";
		}
		filestr<<endl;
	}
	
	filestr.close();
}
double DAG::getPOSPROB(std::string nname){
	double outd=-1.0;
	if(GtoXMap.find(nname)!=GtoXMap.end()){
		int idxlabel = GtoXMap.find(nname)->second;
		if(bpMap.find(nname)!=bpMap.end()){
			factor f = bpMap.find(nname)->second;
			map<int,int > assign;
			assign[idxlabel] =1;
			outd = f(assign);
		}
	}
	return outd;
}

void DAG::bp_revised(string SourceNode,int round){
	
	map<string, factor > MSGFactorMap;
	map<string, factor > NodeFactorMap;
	map<string, factor > EdgeFactorMap;		
			
	
	
	//read Likely hood vector
	vector<double> lhvec;
	fstream filestr;
	string likelyhoodfname ="../data/test2/likelihood";
	filestr.open(likelyhoodfname.c_str(), ios_base::in);
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			lhvec.push_back(atof(str.c_str()));
		}
	}
	
	filestr.close();
	
	 
	map<string,double>::iterator nodeitr;
	int nodeidx =1;	 
	for(nodeitr = nodeMap.begin();nodeitr != nodeMap.end();nodeitr++){
		//convert node factor
		GtoXMap.insert(make_pair(nodeitr->first,nodeidx));
		
		double qvalue = nodeitr->second;
		//check LH
		if (qvalue < 0.000001) {
			qvalue = 0.000001; // fixed truncating error
		}
		
		double likehood= (1-qvalue)/qvalue;
//		cout << likehood << endl;
		if(likehood>1){
			
			map<int,int> scope,assign;
			scope[nodeidx]=2;
			factor nodef(scope);
			assign[nodeidx]=0;
			nodef(assign)=1.0;
			assign[nodeidx]=1;
			nodef(assign)=likehood;
			NodeFactorMap.insert(make_pair(nodeitr->first,nodef));
			//nodef.print(cout);
		}
		else{
			if(likehood == 0){
				likehood = 0.000001; // fixed truncating error
			}
			likehood = 1/likehood;
			map<int,int> scope,assign;
			scope[nodeidx]=2;
			factor nodef(scope);
			assign[nodeidx]=0;
			nodef(assign)=likehood;
			assign[nodeidx]=1;
			nodef(assign)=1.0;
			NodeFactorMap.insert(make_pair(nodeitr->first,nodef));
			//nodef.print(cout);
		}
		
		
				
		nodeidx++;
	}
	
	//convert edge factor
	for(nodeitr = nodeMap.begin();nodeitr != nodeMap.end();nodeitr++){
		string sstr=nodeitr->first;
		map<string,double> LinkMap = edgeMap.find(nodeitr->first)->second;
		map<string,double>::iterator edgeitr;
		for(edgeitr =  LinkMap.begin();edgeitr !=  LinkMap.end();edgeitr++){
			
			double edgeweight = pow(2, edgeitr->second);
			string tstr = edgeitr->first;
			string edgename =  sstr+ ","+ tstr;
			string r_edgename =  tstr+ ","+ sstr;
			map<int,int> edgescope;
			int sidx = GtoXMap.find(sstr)->second;
			int tidx = GtoXMap.find(tstr)->second;
			edgescope[sidx]=2;
			edgescope[tidx]=2;
			factor edgefactor(edgescope);
			map<int,int> assign;
			assign[sidx]=0; assign[tidx]=0;
			edgefactor(assign) = edgeweight;
			assign[sidx]=0; assign[tidx]=1;
			edgefactor(assign) = 1.0;
			assign[sidx]=1; assign[tidx]=0;
			edgefactor(assign) = 1.0;
			assign[sidx]=1; assign[tidx]=1;
			edgefactor(assign) = edgeweight;
			
			if(EdgeFactorMap.find(edgename)==EdgeFactorMap.end()){
				EdgeFactorMap.insert(make_pair(edgename,edgefactor));
				//edgefactor.print(cout);
			}
			if(EdgeFactorMap.find(r_edgename)==EdgeFactorMap.end()){
				EdgeFactorMap.insert(make_pair(r_edgename,edgefactor));
				//edgefactor.print(cout);
				
			}
			//init msgfactor
			map<int,int> msgscope;
			int nodeidx = GtoXMap.find(sstr)->second;
			msgscope[nodeidx]=2;
			factor msgfactor(msgscope);
			map<int,int> msgassin;
			msgassin[nodeidx]=0;
			msgfactor(msgassin) = 1.0;
			msgassin[nodeidx]=1;
			msgfactor(msgassin) = 1.0;
			if(MSGFactorMap.find(r_edgename)==MSGFactorMap.end()){
				MSGFactorMap.insert(make_pair(r_edgename,msgfactor));
				//msgfactor.print(cout);
				
			}
		}
		 
	}
	
	//determine the propagate order DFS postorder
	//get DFS order
	vector<string> DFSVec ; 
	list<string> Stack;
	map<string,string> VisitedMap; 
	Stack.push_front(SourceNode);
	VisitedMap.insert(make_pair(SourceNode,""));
	while(Stack.size()>0){
		string topstr = Stack.front();
		map<string, double> LinkMap= edgeMap.find(topstr)->second;
		string nextstr ="#";
		map<string, double>::iterator edgeitr;
		
		for(edgeitr= LinkMap.begin();edgeitr!=LinkMap.end();edgeitr++){
			if(VisitedMap.find(edgeitr->first)== VisitedMap.end()){
				nextstr = edgeitr->first;
			}
				
		} 
		
		if(nextstr.compare("#")!=0){
			Stack.push_front(nextstr);
			VisitedMap.insert(make_pair(nextstr,""));
		}
		else{
			DFSVec.push_back(topstr);
			Stack.pop_front();
		}
	}
	//
	
	//init oldprob
	map<string, double> oldProbMap;
	
	for(int i =0;i<DFSVec.size();i++){
			string sstr = DFSVec[i];
			oldProbMap.insert(make_pair(sstr,-10));
	}
	
	
	for(int r=0;r<round;r++){
	
	int diff_num =0;
	double diff_v =0.0;
	for(int i =0;i<DFSVec.size();i++){
			string sstr = DFSVec[i];
			double old_prob = oldProbMap.find(sstr)->second;
			if(old_prob != getPOSPROB(sstr)){
			
				diff_num ++;
				diff_v += fabs(old_prob-getPOSPROB(sstr));
			}
			oldProbMap[sstr] = getPOSPROB(sstr);
			
	}	
	double avgdiff = diff_v/(double) diff_num;
	
//  cout << "bp:"<< r << ","<< avgdiff << "," << diff_num<<endl; 
	if(avgdiff<0.001 || diff_num == 0){
		break;
	}
	
	//propagate forward
	
	for(int i =0;i<DFSVec.size();i++){
		//cout << "DFS:"<<DFSVec[i]<<endl;
		string sstr = DFSVec[i];
		map<string,double> LinkMap = edgeMap.find(DFSVec[i])->second;
		//update msg
		map<string,double>::iterator linkitr;
		int sidx =GtoXMap.find(DFSVec[i])->second;	
		map<int,int> msgscope;
		
		msgscope[sidx] =2;
		
		
		for(linkitr = LinkMap.begin();linkitr != LinkMap.end();linkitr++){
			string newmsgname = DFSVec[i]+","+linkitr->first;
			string tstr = linkitr->first;
		
			map<string,double>::iterator initr;
			factor msgfactor(msgscope);	
			map<int,int> msgassign ;
			msgassign[sidx]=0;
			msgfactor(msgassign) = 1.0;
			msgassign[sidx]=1;
			msgfactor(msgassign) = 1.0;
			//multiply incoming msg	
			for(initr = LinkMap.begin();initr!=LinkMap.end();initr++ ){
				string instr = initr->first ;
				
				if(instr.compare(tstr)!=0){
					
					string multistr = instr+","+DFSVec[i];
					
					 msgfactor = msgfactor * MSGFactorMap.find(multistr)->second;
					
				}
			}
			
			//msgfactor.print(cout);
			//
			factor newmsg =  ( EdgeFactorMap.find(newmsgname)->second*NodeFactorMap.find(sstr)->second *msgfactor);
			int tidx = GtoXMap.find( tstr)->second;
			map<int,int> sumoutscope;
			sumoutscope[sidx]=2;
			factor updatedmsg = newmsg.sum(sumoutscope);
			map<int,int> upmsgassign;
			upmsgassign[tidx] =0;
			double nd = 0;
			nd += updatedmsg(upmsgassign);
			upmsgassign[tidx] =1;
			nd += updatedmsg(upmsgassign);
			updatedmsg /= nd;
			
			MSGFactorMap.find(newmsgname)->second = updatedmsg;
			//updatedmsg.print(cout);			
			
		}
		
	}
	//propagate backward
	for(int i =DFSVec.size()-1;i>=0;i--){

		string sstr = DFSVec[i];
		map<string,double> LinkMap = edgeMap.find(DFSVec[i])->second;
		//update msg
		map<string,double>::iterator linkitr;
		int sidx =GtoXMap.find(DFSVec[i])->second;	
		map<int,int> msgscope;
		
		msgscope[sidx] =2;
		
		
		for(linkitr = LinkMap.begin();linkitr != LinkMap.end();linkitr++){
			string newmsgname = DFSVec[i]+","+linkitr->first;
			string tstr = linkitr->first;
			
			map<string,double>::iterator initr;
			factor msgfactor(msgscope);	
			map<int,int> msgassign ;
			msgassign[sidx]=0;
			msgfactor(msgassign) = 1.0;
			msgassign[sidx]=1;
			msgfactor(msgassign) = 1.0;
			//multiply incoming msg	
			for(initr = LinkMap.begin();initr!=LinkMap.end();initr++ ){
				string instr = initr->first ;
				
				if(instr.compare(tstr)!=0){
					
					string multistr = instr+","+DFSVec[i];
					
					 msgfactor = msgfactor * MSGFactorMap.find(multistr)->second;
					
				}
			}
			
			//msgfactor.print(cout);
			//
			factor newmsg =  ( EdgeFactorMap.find(newmsgname)->second*NodeFactorMap.find(sstr)->second *msgfactor);
			int tidx = GtoXMap.find( tstr)->second;
			map<int,int> sumoutscope;
			sumoutscope[sidx]=2;
			factor updatedmsg = newmsg.sum(sumoutscope);
			map<int,int> upmsgassign;
			upmsgassign[tidx] =0;
			double nd = 0;
			nd += updatedmsg(upmsgassign);
			upmsgassign[tidx] =1;
			nd += updatedmsg(upmsgassign);
			updatedmsg /= nd;
			MSGFactorMap.find(newmsgname)->second = updatedmsg;
			//updatedmsg.print(cout);			
			
		}
		//calculate BP
		for(int i =0;i<DFSVec.size();i++){
			string sstr = DFSVec[i];
			map<string,double> LinkMap = edgeMap.find(DFSVec[i])->second;
			//update msg
			map<string,double>::iterator linkitr;
			int sidx =GtoXMap.find(DFSVec[i])->second;	
			map<int,int> msgscope;
			msgscope[sidx] =2;
			factor msgfactor(msgscope);	
			map<int,int> msgassign ;
			msgassign[sidx]=0;
			msgfactor(msgassign) = 1.0;
			msgassign[sidx]=1;
			msgfactor(msgassign) = 1.0;
		
			for(linkitr = LinkMap.begin();linkitr != LinkMap.end();linkitr++){
			
				string newmsgname = DFSVec[i]+","+linkitr->first;
				string tstr = linkitr->first;
		
			//multiply incoming msg	

				string multistr = tstr+","+DFSVec[i];
				msgfactor = msgfactor * MSGFactorMap.find(multistr)->second;
				
			}
			factor updatedmsg = msgfactor*NodeFactorMap.find(sstr)->second;
			map<int,int> upmsgassign;
				upmsgassign[sidx] =0;
				double nd = 0;
				nd += updatedmsg(upmsgassign);
				upmsgassign[sidx] =1;
				nd += updatedmsg(upmsgassign);
				updatedmsg /= nd;
				if(bpMap.find(sstr)==bpMap.end()){
					bpMap.insert(make_pair(sstr, updatedmsg));
				}
				else{
					bpMap.erase(sstr);
					bpMap.insert(make_pair(sstr, updatedmsg));
				}
			//cout<< "bp "<< sstr<<endl;
			//updatedmsg.print(cout);
		}
	}
	
	
	}
	
	
	
	
}

void DAG::bp_noiseq(string SourceNode,int round){
	
	map<string, factor > MSGFactorMap;
	map<string, factor > NodeFactorMap;
	map<string, factor > EdgeFactorMap;		
			
	
	
	//read Likely hood vector
	vector<double> lhvec;
	fstream filestr;
	string likelyhoodfname ="../data/test2/likelihood";
	filestr.open(likelyhoodfname.c_str(), ios_base::in);
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			lhvec.push_back(atof(str.c_str()));
		}
	}
	
	filestr.close();
	
	 
	map<string,double>::iterator nodeitr;
	int nodeidx =1;	 
	for(nodeitr = nodeMap.begin();nodeitr != nodeMap.end();nodeitr++){
		//convert node factor
		GtoXMap.insert(make_pair(nodeitr->first,nodeidx));
		
		double qvalue = nodeitr->second;
		//check LH
		if (qvalue < 0.000001) {
			qvalue = 0.000001; // fixed truncating error
		}
		double likehood;
		if((1-qvalue) ==0){
			likehood = 1000000;
		}
		else{
			likehood = (qvalue)/(1-qvalue);
		}
//		cout << likehood << endl;
		if(likehood>1){
			
			map<int,int> scope,assign;
			scope[nodeidx]=2;
			factor nodef(scope);
			assign[nodeidx]=0;
			nodef(assign)=1.0;
			assign[nodeidx]=1;
			nodef(assign)=likehood;
			NodeFactorMap.insert(make_pair(nodeitr->first,nodef));
			//nodef.print(cout);
		}
		else{
			if(likehood == 0){
				likehood = 0.000001; // fixed truncating error
			}
			likehood = 1/likehood;
			map<int,int> scope,assign;
			scope[nodeidx]=2;
			factor nodef(scope);
			assign[nodeidx]=0;
			nodef(assign)=likehood;
			assign[nodeidx]=1;
			nodef(assign)=1.0;
			NodeFactorMap.insert(make_pair(nodeitr->first,nodef));
			//nodef.print(cout);
		}
		
		
				
		nodeidx++;
	}
	
	//convert edge factor
	for(nodeitr = nodeMap.begin();nodeitr != nodeMap.end();nodeitr++){
		string sstr=nodeitr->first;
		map<string,double> LinkMap = edgeMap.find(nodeitr->first)->second;
		map<string,double>::iterator edgeitr;
		for(edgeitr =  LinkMap.begin();edgeitr !=  LinkMap.end();edgeitr++){
			
			double edgeweight = pow(2, edgeitr->second);
			string tstr = edgeitr->first;
			string edgename =  sstr+ ","+ tstr;
			string r_edgename =  tstr+ ","+ sstr;
			map<int,int> edgescope;
			int sidx = GtoXMap.find(sstr)->second;
			int tidx = GtoXMap.find(tstr)->second;
			edgescope[sidx]=2;
			edgescope[tidx]=2;
			factor edgefactor(edgescope);
			map<int,int> assign;
			assign[sidx]=0; assign[tidx]=0;
			edgefactor(assign) = edgeweight;
			assign[sidx]=0; assign[tidx]=1;
			edgefactor(assign) = 1.0;
			assign[sidx]=1; assign[tidx]=0;
			edgefactor(assign) = 1.0;
			assign[sidx]=1; assign[tidx]=1;
			edgefactor(assign) = edgeweight;
			
			if(EdgeFactorMap.find(edgename)==EdgeFactorMap.end()){
				EdgeFactorMap.insert(make_pair(edgename,edgefactor));
				//edgefactor.print(cout);
			}
			if(EdgeFactorMap.find(r_edgename)==EdgeFactorMap.end()){
				EdgeFactorMap.insert(make_pair(r_edgename,edgefactor));
				//edgefactor.print(cout);
				
			}
			//init msgfactor
			map<int,int> msgscope;
			int nodeidx = GtoXMap.find(sstr)->second;
			msgscope[nodeidx]=2;
			factor msgfactor(msgscope);
			map<int,int> msgassin;
			msgassin[nodeidx]=0;
			msgfactor(msgassin) = 1.0;
			msgassin[nodeidx]=1;
			msgfactor(msgassin) = 1.0;
			if(MSGFactorMap.find(r_edgename)==MSGFactorMap.end()){
				MSGFactorMap.insert(make_pair(r_edgename,msgfactor));
				//msgfactor.print(cout);
				
			}
		}
		 
	}
	
	//determine the propagate order DFS postorder
	//get DFS order
	vector<string> DFSVec ; 
	list<string> Stack;
	map<string,string> VisitedMap; 
	Stack.push_front(SourceNode);
	VisitedMap.insert(make_pair(SourceNode,""));
	while(Stack.size()>0){
		string topstr = Stack.front();
		map<string, double> LinkMap= edgeMap.find(topstr)->second;
		string nextstr ="#";
		map<string, double>::iterator edgeitr;
		
		for(edgeitr= LinkMap.begin();edgeitr!=LinkMap.end();edgeitr++){
			if(VisitedMap.find(edgeitr->first)== VisitedMap.end()){
				nextstr = edgeitr->first;
			}
				
		} 
		
		if(nextstr.compare("#")!=0){
			Stack.push_front(nextstr);
			VisitedMap.insert(make_pair(nextstr,""));
		}
		else{
			DFSVec.push_back(topstr);
			Stack.pop_front();
		}
	}
	//
	
	//init oldprob
	map<string, double> oldProbMap;
	
	for(int i =0;i<DFSVec.size();i++){
			string sstr = DFSVec[i];
			oldProbMap.insert(make_pair(sstr,-10));
	}
	
	
	for(int r=0;r<round;r++){
	
	int diff_num =0;
	double diff_v =0.0;
	for(int i =0;i<DFSVec.size();i++){
			string sstr = DFSVec[i];
			double old_prob = oldProbMap.find(sstr)->second;
			if(old_prob != getPOSPROB(sstr)){
			
				diff_num ++;
				diff_v += fabs(old_prob-getPOSPROB(sstr));
			}
			oldProbMap[sstr] = getPOSPROB(sstr);
			
	}	
	double avgdiff = diff_v/(double) diff_num;
	
//  cout << "bp:"<< r << ","<< avgdiff << "," << diff_num<<endl; 
	if(avgdiff<0.001 || diff_num == 0){
		break;
	}
	
	//propagate forward
	
	for(int i =0;i<DFSVec.size();i++){
		//cout << "DFS:"<<DFSVec[i]<<endl;
		string sstr = DFSVec[i];
		map<string,double> LinkMap = edgeMap.find(DFSVec[i])->second;
		//update msg
		map<string,double>::iterator linkitr;
		int sidx =GtoXMap.find(DFSVec[i])->second;	
		map<int,int> msgscope;
		
		msgscope[sidx] =2;
		
		
		for(linkitr = LinkMap.begin();linkitr != LinkMap.end();linkitr++){
			string newmsgname = DFSVec[i]+","+linkitr->first;
			string tstr = linkitr->first;
		
			map<string,double>::iterator initr;
			factor msgfactor(msgscope);	
			map<int,int> msgassign ;
			msgassign[sidx]=0;
			msgfactor(msgassign) = 1.0;
			msgassign[sidx]=1;
			msgfactor(msgassign) = 1.0;
			//multiply incoming msg	
			for(initr = LinkMap.begin();initr!=LinkMap.end();initr++ ){
				string instr = initr->first ;
				
				if(instr.compare(tstr)!=0){
					
					string multistr = instr+","+DFSVec[i];
					
					 msgfactor = msgfactor * MSGFactorMap.find(multistr)->second;
					
				}
			}
			
			//msgfactor.print(cout);
			//
			factor newmsg =  ( EdgeFactorMap.find(newmsgname)->second*NodeFactorMap.find(sstr)->second *msgfactor);
			int tidx = GtoXMap.find( tstr)->second;
			map<int,int> sumoutscope;
			sumoutscope[sidx]=2;
			factor updatedmsg = newmsg.sum(sumoutscope);
			map<int,int> upmsgassign;
			upmsgassign[tidx] =0;
			double nd = 0;
			nd += updatedmsg(upmsgassign);
			upmsgassign[tidx] =1;
			nd += updatedmsg(upmsgassign);
			updatedmsg /= nd;
			
			MSGFactorMap.find(newmsgname)->second = updatedmsg;
			//updatedmsg.print(cout);			
			
		}
		
	}
	//propagate backward
	for(int i =DFSVec.size()-1;i>=0;i--){

		string sstr = DFSVec[i];
		map<string,double> LinkMap = edgeMap.find(DFSVec[i])->second;
		//update msg
		map<string,double>::iterator linkitr;
		int sidx =GtoXMap.find(DFSVec[i])->second;	
		map<int,int> msgscope;
		
		msgscope[sidx] =2;
		
		
		for(linkitr = LinkMap.begin();linkitr != LinkMap.end();linkitr++){
			string newmsgname = DFSVec[i]+","+linkitr->first;
			string tstr = linkitr->first;
			
			map<string,double>::iterator initr;
			factor msgfactor(msgscope);	
			map<int,int> msgassign ;
			msgassign[sidx]=0;
			msgfactor(msgassign) = 1.0;
			msgassign[sidx]=1;
			msgfactor(msgassign) = 1.0;
			//multiply incoming msg	
			for(initr = LinkMap.begin();initr!=LinkMap.end();initr++ ){
				string instr = initr->first ;
				
				if(instr.compare(tstr)!=0){
					
					string multistr = instr+","+DFSVec[i];
					
					 msgfactor = msgfactor * MSGFactorMap.find(multistr)->second;
					
				}
			}
			
			//msgfactor.print(cout);
			//
			factor newmsg =  ( EdgeFactorMap.find(newmsgname)->second*NodeFactorMap.find(sstr)->second *msgfactor);
			int tidx = GtoXMap.find( tstr)->second;
			map<int,int> sumoutscope;
			sumoutscope[sidx]=2;
			factor updatedmsg = newmsg.sum(sumoutscope);
			map<int,int> upmsgassign;
			upmsgassign[tidx] =0;
			double nd = 0;
			nd += updatedmsg(upmsgassign);
			upmsgassign[tidx] =1;
			nd += updatedmsg(upmsgassign);
			updatedmsg /= nd;
			MSGFactorMap.find(newmsgname)->second = updatedmsg;
			//updatedmsg.print(cout);			
			
		}
		//calculate BP
		for(int i =0;i<DFSVec.size();i++){
			string sstr = DFSVec[i];
			map<string,double> LinkMap = edgeMap.find(DFSVec[i])->second;
			//update msg
			map<string,double>::iterator linkitr;
			int sidx =GtoXMap.find(DFSVec[i])->second;	
			map<int,int> msgscope;
			msgscope[sidx] =2;
			factor msgfactor(msgscope);	
			map<int,int> msgassign ;
			msgassign[sidx]=0;
			msgfactor(msgassign) = 1.0;
			msgassign[sidx]=1;
			msgfactor(msgassign) = 1.0;
		
			for(linkitr = LinkMap.begin();linkitr != LinkMap.end();linkitr++){
			
				string newmsgname = DFSVec[i]+","+linkitr->first;
				string tstr = linkitr->first;
		
			//multiply incoming msg	

				string multistr = tstr+","+DFSVec[i];
				msgfactor = msgfactor * MSGFactorMap.find(multistr)->second;
				
			}
			factor updatedmsg = msgfactor*NodeFactorMap.find(sstr)->second;
			map<int,int> upmsgassign;
				upmsgassign[sidx] =0;
				double nd = 0;
				nd += updatedmsg(upmsgassign);
				upmsgassign[sidx] =1;
				nd += updatedmsg(upmsgassign);
				updatedmsg /= nd;
				if(bpMap.find(sstr)==bpMap.end()){
					bpMap.insert(make_pair(sstr, updatedmsg));
				}
				else{
					bpMap.erase(sstr);
					bpMap.insert(make_pair(sstr, updatedmsg));
				}
			//cout<< "bp "<< sstr<<endl;
			//updatedmsg.print(cout);
		}
	}
	
	
	}
	
	
	
	
}

void DAG::bp(string SourceNode,int round){
	
	map<string, factor > MSGFactorMap;
	map<string, factor > NodeFactorMap;
	map<string, factor > EdgeFactorMap;		
			
	
	
	//read Likely hood vector
	vector<double> lhvec;
	fstream filestr;
	string likelyhoodfname ="../data/test2/likelihood";
	filestr.open(likelyhoodfname.c_str(), ios_base::in);
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			lhvec.push_back(atof(str.c_str()));
		}
	}
	
	filestr.close();
	
	 
	map<string,double>::iterator nodeitr;
	int nodeidx =1;	 
	for(nodeitr = nodeMap.begin();nodeitr != nodeMap.end();nodeitr++){
		//convert node factor
		GtoXMap.insert(make_pair(nodeitr->first,nodeidx));
		
		double qvalue = nodeitr->second;
		//check LH
		int idx = (qvalue*100)/5;
		if (idx ==20){
			idx =19;
		}
		double likehood;
		if(lhvec[idx]>1){
			likehood = lhvec[idx];
			map<int,int> scope,assign;
			scope[nodeidx]=2;
			factor nodef(scope);
			assign[nodeidx]=0;
			nodef(assign)=1.0;
			assign[nodeidx]=1;
			nodef(assign)=likehood;
			NodeFactorMap.insert(make_pair(nodeitr->first,nodef));
			//nodef.print(cout);
		}
		else{
			likehood = 1/(double)lhvec[idx];
			map<int,int> scope,assign;
			scope[nodeidx]=2;
			factor nodef(scope);
			assign[nodeidx]=0;
			nodef(assign)=likehood;
			assign[nodeidx]=1;
			nodef(assign)=1.0;
			NodeFactorMap.insert(make_pair(nodeitr->first,nodef));
			//nodef.print(cout);
		}
		
		
				
		nodeidx++;
	}
	
	//convert edge factor
	for(nodeitr = nodeMap.begin();nodeitr != nodeMap.end();nodeitr++){
		string sstr=nodeitr->first;
		map<string,double> LinkMap = edgeMap.find(nodeitr->first)->second;
		map<string,double>::iterator edgeitr;
		for(edgeitr =  LinkMap.begin();edgeitr !=  LinkMap.end();edgeitr++){
			
			double edgeweight = pow(2, edgeitr->second);
			string tstr = edgeitr->first;
			string edgename =  sstr+ ","+ tstr;
			string r_edgename =  tstr+ ","+ sstr;
			map<int,int> edgescope;
			int sidx = GtoXMap.find(sstr)->second;
			int tidx = GtoXMap.find(tstr)->second;
			edgescope[sidx]=2;
			edgescope[tidx]=2;
			factor edgefactor(edgescope);
			map<int,int> assign;
			assign[sidx]=0; assign[tidx]=0;
			edgefactor(assign) = edgeweight;
			assign[sidx]=0; assign[tidx]=1;
			edgefactor(assign) = 1.0;
			assign[sidx]=1; assign[tidx]=0;
			edgefactor(assign) = 1.0;
			assign[sidx]=1; assign[tidx]=1;
			edgefactor(assign) = edgeweight;
			
			if(EdgeFactorMap.find(edgename)==EdgeFactorMap.end()){
				EdgeFactorMap.insert(make_pair(edgename,edgefactor));
				//edgefactor.print(cout);
			}
			if(EdgeFactorMap.find(r_edgename)==EdgeFactorMap.end()){
				EdgeFactorMap.insert(make_pair(r_edgename,edgefactor));
				//edgefactor.print(cout);
				
			}
			//init msgfactor
			map<int,int> msgscope;
			int nodeidx = GtoXMap.find(sstr)->second;
			msgscope[nodeidx]=2;
			factor msgfactor(msgscope);
			map<int,int> msgassin;
			msgassin[nodeidx]=0;
			msgfactor(msgassin) = 1.0;
			msgassin[nodeidx]=1;
			msgfactor(msgassin) = 1.0;
			if(MSGFactorMap.find(r_edgename)==MSGFactorMap.end()){
				MSGFactorMap.insert(make_pair(r_edgename,msgfactor));
				//msgfactor.print(cout);
				
			}
		}
		 
	}
	
	//determine the propagate order DFS postorder
	//get DFS order
	vector<string> DFSVec ; 
	list<string> Stack;
	map<string,string> VisitedMap; 
	Stack.push_front(SourceNode);
	VisitedMap.insert(make_pair(SourceNode,""));
	while(Stack.size()>0){
		string topstr = Stack.front();
		map<string, double> LinkMap= edgeMap.find(topstr)->second;
		string nextstr ="#";
		map<string, double>::iterator edgeitr;
		
		for(edgeitr= LinkMap.begin();edgeitr!=LinkMap.end();edgeitr++){
			if(VisitedMap.find(edgeitr->first)== VisitedMap.end()){
				nextstr = edgeitr->first;
			}
				
		} 
		
		if(nextstr.compare("#")!=0){
			Stack.push_front(nextstr);
			VisitedMap.insert(make_pair(nextstr,""));
		}
		else{
			DFSVec.push_back(topstr);
			Stack.pop_front();
		}
	}
	//
	
	//init oldprob
	map<string, double> oldProbMap;
	
	for(int i =0;i<DFSVec.size();i++){
			string sstr = DFSVec[i];
			oldProbMap.insert(make_pair(sstr,-10));
	}
	
	
	for(int r=0;r<round;r++){
	
	int diff_num =0;
	double diff_v =0.0;
	for(int i =0;i<DFSVec.size();i++){
			string sstr = DFSVec[i];
			double old_prob = oldProbMap.find(sstr)->second;
			if(old_prob != getPOSPROB(sstr)){
			
				diff_num ++;
				diff_v += fabs(old_prob-getPOSPROB(sstr));
			}
			oldProbMap[sstr] = getPOSPROB(sstr);
			
	}	
	double avgdiff = diff_v/(double) diff_num;
	
  //cout << "bp:"<< r << ","<< avgdiff << "," << diff_num<<endl; 
	if(avgdiff<0.001 || diff_num == 0){
		break;
	}
	
	//propagate forward
	
	for(int i =0;i<DFSVec.size();i++){
		//cout << "DFS:"<<DFSVec[i]<<endl;
		string sstr = DFSVec[i];
		map<string,double> LinkMap = edgeMap.find(DFSVec[i])->second;
		//update msg
		map<string,double>::iterator linkitr;
		int sidx =GtoXMap.find(DFSVec[i])->second;	
		map<int,int> msgscope;
		
		msgscope[sidx] =2;
		
		
		for(linkitr = LinkMap.begin();linkitr != LinkMap.end();linkitr++){
			string newmsgname = DFSVec[i]+","+linkitr->first;
			string tstr = linkitr->first;
		
			map<string,double>::iterator initr;
			factor msgfactor(msgscope);	
			map<int,int> msgassign ;
			msgassign[sidx]=0;
			msgfactor(msgassign) = 1.0;
			msgassign[sidx]=1;
			msgfactor(msgassign) = 1.0;
			//multiply incoming msg	
			for(initr = LinkMap.begin();initr!=LinkMap.end();initr++ ){
				string instr = initr->first ;
				
				if(instr.compare(tstr)!=0){
					
					string multistr = instr+","+DFSVec[i];
					
					 msgfactor = msgfactor * MSGFactorMap.find(multistr)->second;
					
				}
			}
			
			//msgfactor.print(cout);
			//
			factor newmsg =  ( EdgeFactorMap.find(newmsgname)->second*NodeFactorMap.find(sstr)->second *msgfactor);
			int tidx = GtoXMap.find( tstr)->second;
			map<int,int> sumoutscope;
			sumoutscope[sidx]=2;
			factor updatedmsg = newmsg.sum(sumoutscope);
			map<int,int> upmsgassign;
			upmsgassign[tidx] =0;
			double nd = 0;
			nd += updatedmsg(upmsgassign);
			upmsgassign[tidx] =1;
			nd += updatedmsg(upmsgassign);
			updatedmsg /= nd;
			
			MSGFactorMap.find(newmsgname)->second = updatedmsg;
			//updatedmsg.print(cout);			
			
		}
		
	}
	//propagate backward
	for(int i =DFSVec.size()-1;i>=0;i--){

		string sstr = DFSVec[i];
		map<string,double> LinkMap = edgeMap.find(DFSVec[i])->second;
		//update msg
		map<string,double>::iterator linkitr;
		int sidx =GtoXMap.find(DFSVec[i])->second;	
		map<int,int> msgscope;
		
		msgscope[sidx] =2;
		
		
		for(linkitr = LinkMap.begin();linkitr != LinkMap.end();linkitr++){
			string newmsgname = DFSVec[i]+","+linkitr->first;
			string tstr = linkitr->first;
			
			map<string,double>::iterator initr;
			factor msgfactor(msgscope);	
			map<int,int> msgassign ;
			msgassign[sidx]=0;
			msgfactor(msgassign) = 1.0;
			msgassign[sidx]=1;
			msgfactor(msgassign) = 1.0;
			//multiply incoming msg	
			for(initr = LinkMap.begin();initr!=LinkMap.end();initr++ ){
				string instr = initr->first ;
				
				if(instr.compare(tstr)!=0){
					
					string multistr = instr+","+DFSVec[i];
					
					 msgfactor = msgfactor * MSGFactorMap.find(multistr)->second;
					
				}
			}
			
			//msgfactor.print(cout);
			//
			factor newmsg =  ( EdgeFactorMap.find(newmsgname)->second*NodeFactorMap.find(sstr)->second *msgfactor);
			int tidx = GtoXMap.find( tstr)->second;
			map<int,int> sumoutscope;
			sumoutscope[sidx]=2;
			factor updatedmsg = newmsg.sum(sumoutscope);
			map<int,int> upmsgassign;
			upmsgassign[tidx] =0;
			double nd = 0;
			nd += updatedmsg(upmsgassign);
			upmsgassign[tidx] =1;
			nd += updatedmsg(upmsgassign);
			updatedmsg /= nd;
			MSGFactorMap.find(newmsgname)->second = updatedmsg;
			//updatedmsg.print(cout);			
			
		}
		//calculate BP
		for(int i =0;i<DFSVec.size();i++){
			string sstr = DFSVec[i];
			map<string,double> LinkMap = edgeMap.find(DFSVec[i])->second;
			//update msg
			map<string,double>::iterator linkitr;
			int sidx =GtoXMap.find(DFSVec[i])->second;	
			map<int,int> msgscope;
			msgscope[sidx] =2;
			factor msgfactor(msgscope);	
			map<int,int> msgassign ;
			msgassign[sidx]=0;
			msgfactor(msgassign) = 1.0;
			msgassign[sidx]=1;
			msgfactor(msgassign) = 1.0;
		
			for(linkitr = LinkMap.begin();linkitr != LinkMap.end();linkitr++){
			
				string newmsgname = DFSVec[i]+","+linkitr->first;
				string tstr = linkitr->first;
		
			//multiply incoming msg	

				string multistr = tstr+","+DFSVec[i];
				msgfactor = msgfactor * MSGFactorMap.find(multistr)->second;
				
			}
			factor updatedmsg = msgfactor*NodeFactorMap.find(sstr)->second;
			map<int,int> upmsgassign;
				upmsgassign[sidx] =0;
				double nd = 0;
				nd += updatedmsg(upmsgassign);
				upmsgassign[sidx] =1;
				nd += updatedmsg(upmsgassign);
				updatedmsg /= nd;
				if(bpMap.find(sstr)==bpMap.end()){
					bpMap.insert(make_pair(sstr, updatedmsg));
				}
				else{
					bpMap.erase(sstr);
					bpMap.insert(make_pair(sstr, updatedmsg));
				}
			//cout<< "bp "<< sstr<<endl;
			//updatedmsg.print(cout);
		}
	}
	
	
	}
	
	
	
	
}



//vector<string> DAG::scheduleNodes(string SourceNode){
//	  vector<string> OutVec;
//	
//		list<string > Q;		
//		map<string,double> addednodesMap;
//		
//		Q.push_back(SourceNode);
//		addednodesMap.insert(make_pair(SourceNode,0.0));
//		map<string, double>::iterator edge_itr;
//		while(Q.size()>0){
//			
//			string sstr = Q.front();
//	 		Q.pop_front();
//	 		
//	 		//cout<< "add " <<sstr <<endl;
//	 	 //sstr is the current node.
//		 	edgeitr = edgeMap.find(sstr);
//		 	map<string,double> neighborMap = edgeitr->second;
//		 
//		 	for(neighboritr = neighborMap.begin();neighboritr !=neighborMap.end();neighboritr++){
//		 		addeditr = addednodesMap.find(neighboritr->first);
//		 		if(addeditr==addednodesMap.end()){
//		 			
//		 			Q.push_back(neighboritr->first);
//		 			addednodesMap.insert(make_pair(neighboritr->first,0.0));
//		 		}
//		 }
//		 
//		}
//	
//	
//	
//	return OutVec;
//}
void DAG::read(string FName){
	fstream filestr;
	filestr.open(FName.c_str(), ios_base::in);
	string str;	 
	vector<string> tokvec;
	TOK tok;
	while(!filestr.eof()){
		getline(filestr,str);
	 	string NodeStr,EdgeStr;
	 	if(str.length()>1){
	 		vector<string> TokVec = tok.Tokenize(str,":");
			
			
			
	 		NodeStr = TokVec[0];
			
			if(TokVec.size()>1){
	 			EdgeStr = TokVec[1];
			}
			else{
				EdgeStr="";
			}

	 		TokVec = tok.Tokenize(NodeStr,"_");
	 		string NodeName = TokVec[0];
	 		nodeMap.insert(make_pair( NodeName, atof(TokVec[1].c_str()) ));
	 		map<string,double> LinkMap;
			
			if(EdgeStr.length()==0){
				edgeMap.insert(make_pair( NodeName, LinkMap));
			}
			else{

	 			TokVec = tok.Tokenize(EdgeStr,",");
				for(int i =0;i<TokVec.size();i++){
					vector<string> TmpVec = tok.Tokenize(TokVec[i],"_");
					LinkMap.insert(make_pair( TmpVec[0], atof(TmpVec[1].c_str()) )); 
				}
				edgeMap.insert(make_pair( NodeName, LinkMap));
			}		
	 	}
	
	}
//	 	
	filestr.close();
}


DAG DAG::merge(DAG ADAG,DAG BDAG){
	DAG OutDAG;
	
	//merge A
	map<string,double>::iterator nodeitr,findnodeitr ;

	map<string, map<string, double> >::iterator Outedgeitr,myedgeitr;
	
	for(nodeitr = ADAG.nodeMap.begin();nodeitr != ADAG.nodeMap.end(); nodeitr++){
		OutDAG.nodeMap.insert(make_pair(nodeitr->first,nodeitr->second));
		//initail edges from this node
		myedgeitr = ADAG.edgeMap.find(nodeitr->first);
		map<string,double> LinkMap= myedgeitr->second;
		OutDAG.edgeMap.insert( make_pair(nodeitr->first, LinkMap) );
		
		
	}
	
	for(nodeitr = BDAG.nodeMap.begin();nodeitr != BDAG.nodeMap.end(); nodeitr++){
		findnodeitr = OutDAG.nodeMap.find(nodeitr->first);
		if(findnodeitr == OutDAG.nodeMap.end()){ // new to OutDag
			OutDAG.nodeMap.insert(make_pair(nodeitr->first,nodeitr->second));
			//initail edges from this node
			myedgeitr = BDAG.edgeMap.find(nodeitr->first);
			map<string,double> LinkMap= myedgeitr->second;
			OutDAG.edgeMap.insert( make_pair(nodeitr->first, LinkMap) );

		}
		else{ //OutDag contain  
			Outedgeitr = OutDAG.edgeMap.find(nodeitr->first);
			
			myedgeitr = BDAG.edgeMap.find(nodeitr->first);
			map<string,double> BLinkMap= myedgeitr->second;
			map<string,double>::iterator BLinkItr ;
				for(BLinkItr  = BLinkMap.begin();BLinkItr  != BLinkMap.end();BLinkItr++){
					if(Outedgeitr->second.find(BLinkItr->first) == Outedgeitr->second.end()){
						Outedgeitr->second.insert(make_pair(BLinkItr->first,BLinkItr->second) );
					}
				}
			
		}
		
	}
	
	
	
//	for(edgeitr = ADAG.edgeMap.begin();edgeitr != ADAG.edgeMap.end(); edgeitr++){
//		OutDAG.edgeMap.insert(make_pair(edgeitr->first,edgeitr->second));
//	}
//	
//	for(edgeitr = BDAG.edgeMap.begin();edgeitr != BDAG.edgeMap.end(); edgeitr++){
//		
//		findedgeitr = OutDAG.edgeMap.find(edgeitr->first);
//		
//		map<string,double> BLinkMap;
//		map<string,double>::iterator BLinkMapItr,OutLinkItr;
//		
//			
//		if(findedgeitr == OutDAG.edgeMap.end()){
//			OutDAG.edgeMap.insert(make_pair(edgeitr->first,edgeitr->second));
//		}
//		else{
//			BLinkMap = edgeitr->second;				
//			for(BLinkMapItr = BLinkMap.begin();BLinkMapItr != BLinkMap.begin();BLinkMapItr++){
//				
//				OutLinkItr = findedgeitr->second.find(BLinkMapItr->first);
//				if(OutLinkItr == findedgeitr->second.end()){
//					findedgeitr->second.insert(make_pair( BLinkMapItr->first,BLinkMapItr->second ));
//				}
//			}
//			
//		}
//	}
//	
	
	return OutDAG;
}
vector<string> DAG::partition(string source, DAG & newdag){
	vector<string> outVec;
	list<string> Q;
	map<string,string> parentMap;
	parentMap.insert(make_pair(source,"*"));
	
	Q.push_back(source);
	map<string, double>::iterator edge_itr;
	while(Q.size()>0){
		string vstr = Q.front();
		Q.pop_front();
		
		if(vstr[0]!='d'){
			outVec.push_back(vstr);
		}
		
		map<string,double> emap= (newdag.edgeMap).find(vstr)->second;
		
		for(edge_itr= emap.begin();edge_itr!=emap.end();edge_itr++){
			string dstr = edge_itr->first;
			double eweight = edge_itr->second;
			
			if((parentMap.find(dstr)==parentMap.end()) && eweight>0 ){
				Q.push_back(dstr);
				parentMap.insert(make_pair(dstr,vstr));
				
			}
		}
		
	}
	return outVec;
}

vector<string> DAG::mincut(string source, string sink){
	vector<string> outVec;
	DAG newdag = antiparaedgefreeDAG().residueG();
	 
	//newdag.print(); 
	vector<string> pathVec= BFS(source,sink,newdag);
  
  //while path not 0
	//find min of p
	
	while(pathVec.size()>1){
		double minw=DBL_MAX;
//  
		for(int i = 0 ;i<pathVec.size()-1;i++){
			string vstr = pathVec[i];
			string dstr = pathVec[i+1];
//		
			map<string, double> emap = newdag.edgeMap.find(vstr)->second;
			double ew = emap.find(dstr)->second;
			if(ew<minw){
				minw=ew;
			}	
	
		}
	
	//augement flow
		for(int i = 0 ;i<pathVec.size()-1;i++){
			string vstr = pathVec[i];
			string dstr = pathVec[i+1];
//		
			map<string, double> & emap = newdag.edgeMap.find(vstr)->second;
			double ew = emap.find(dstr)->second;
			ew = ew-minw;
			emap[dstr]=ew;
		
			map<string, double> & fmap = newdag.edgeMap.find(dstr)->second;
			ew = fmap.find(vstr)->second;
			ew = ew+minw;
			fmap[vstr]=ew;
			
		}	
		
	
		pathVec= BFS(source,sink,newdag);
	}
	//newdag.print(); 
  vector<string> tmppar = partition(source, newdag);
  
 // cout<< "we have " << tmppar.size()-1 << " positive prediction and they are " <<endl;
  
  for(int i =0 ; i<tmppar.size();i++){
  	if( tmppar[i].compare(source)!= 0)
  	{
  		outVec.push_back(tmppar[i]);
  	}
  }

//	cout<< "we have " << nodeMap.size()-2-tmppar.size()+1 << " negative prediction and they are " <<endl;
//	
//	map<string, double>::iterator nitr;
//	for(nitr=nodeMap.begin();nitr!=nodeMap.end();nitr++){
//		int flag =0;
//		
//		for(int i =0 ; i<tmppar.size();i++){
//  		if(nitr->first.compare(tmppar[i])==0){
//  			flag=1;
//  		}
//  	}
//		
//		if(flag==0){
//			cout<<nitr->first<<endl;
//		}
//	}
	
	return outVec;
	
}
void DAG::print(){
		map<string, double>::iterator node_itr,edge_itr;
			
		for(node_itr = this-> nodeMap.begin();node_itr!=this-> nodeMap.end();node_itr++){
			cout<< "node : "<<node_itr->first << " , "<< node_itr->second <<endl; 
			if(this-> edgeMap.find(node_itr->first)!=this-> edgeMap.end()){
				map<string,double> emap = (this-> edgeMap.find(node_itr->first))->second;
				for(edge_itr=emap.begin();edge_itr!=emap.end();edge_itr++ ){
					cout << "edge : (" << node_itr->first << "," << edge_itr->first <<") , " <<edge_itr->second <<endl;
				}
			}
		}	
			
}

vector<string> DAG::BFS(string sstr,string tstr, DAG &inputG){
	vector<string> outVec,tmpVec;
	list<string> Q;
	map<string,string> parentMap;
	parentMap.insert(make_pair(sstr,"*"));
	
	Q.push_back(sstr);
	map<string, double>::iterator edge_itr;
	while(Q.size()>0){
		string vstr = Q.front();
		Q.pop_front();
		map<string,double> emap= (inputG.edgeMap).find(vstr)->second;
		
		for(edge_itr= emap.begin();edge_itr!=emap.end();edge_itr++){
			string dstr = edge_itr->first;
			double eweight = edge_itr->second;
			
			if((parentMap.find(dstr)==parentMap.end()) && eweight>0 ){
				Q.push_back(dstr);
				parentMap.insert(make_pair(dstr,vstr));
				//cout<<vstr << " -> "<< dstr <<endl;
			}
		}
		
	}
	//if path exists....	
	if(parentMap.find(tstr)!=parentMap.end()){
			tmpVec.push_back(tstr);
			string parentstr = tstr;
			 
			while(parentMap.find(parentstr)!=parentMap.end()&&parentMap.find(parentstr)->second!="*"){
				parentstr = parentMap.find(parentstr)->second;
				tmpVec.push_back(parentstr);
								
			}
			
	}
		
		for(int i =tmpVec.size()-1;i>=0;i--){
			outVec.push_back(tmpVec[i]);
			//cout <<"path :" << tmpVec[i] <<endl;
		}
	return outVec;
	
}

DAG DAG::residueG(){
	DAG outDAG;
	outDAG.nodeMap = this->nodeMap;
  //outDAG.edgeMap = this->edgeMap;
  map<string, double>::iterator node_itr,edge_itr;
  
  for(node_itr = outDAG.nodeMap.begin();node_itr!= outDAG.nodeMap.end();node_itr++){
  	map<string,double> emap;
  	outDAG.edgeMap.insert(make_pair(node_itr->first, emap));  	
  }
  
  for(node_itr = this->nodeMap.begin();node_itr!= this->nodeMap.end();node_itr++){
  		if(this->edgeMap.find(node_itr->first)!=this->edgeMap.end()){
  			map<string,double> emap= this->edgeMap.find(node_itr->first)->second;  			
  			for(edge_itr = emap.begin();edge_itr!=emap.end();edge_itr++){
  				(outDAG.edgeMap.find(node_itr->first)->second).insert(make_pair(edge_itr->first,edge_itr->second));
  				(outDAG.edgeMap.find(edge_itr->first)->second).insert(make_pair(node_itr->first,0));
  			}
  		}
  }
 
	return outDAG;
}

DAG DAG::antiparaedgefreeDAG(){
	DAG outDAG;
	stringstream ss ;
	map<string, string> antimap;
	int dummy_num =0;
	//outDAG.nodeMap = this->nodeMap;
	//outDAG.edgeMap = this->edgeMap;
	
	map<string, double>::iterator node_itr,edge_itr;
	
	for(node_itr = this->nodeMap.begin();node_itr!=this->nodeMap.end();node_itr++){
		//cout << node_itr->first << ","<< node_itr->second <<endl;
		string node_id = node_itr->first;
		
		outDAG.nodeMap.insert(make_pair(node_id,node_itr->second));
		map<string,double> newemap;
		if(this->edgeMap.find(node_id)!=this->edgeMap.end()){
			map<string, double> eMap = this->edgeMap.find(node_id)->second;
			for(edge_itr=eMap.begin();edge_itr!=eMap.end();edge_itr++){
				//cout << edge_itr->first << ","<< edge_itr->second <<endl;
				
				string dstr = edge_itr->first; 
				double eweight = edge_itr->second;
				int is_anti =0;
				//detect if this edge is an anti-parallel edge
				if(this->edgeMap.find(dstr)!=this->edgeMap.end()){
					map<string, double> dMap =this->edgeMap.find(dstr)->second;
					ss.str("");
					ss<< dstr<< ","<< node_id;
					if(dMap.find(node_id)!=dMap.end()&& antimap.find(ss.str())==antimap.end()){
						//cout<< " anti para edge "<< node_id << ","<< dstr<< endl;
						ss.str("");
						ss<< node_id << ","<< dstr;;
						antimap.insert(make_pair(ss.str(),""));
						is_anti=1;
					}				
				}
				
				//para edge
				if(is_anti==1){
					ss.str("");
					ss<< "d"<< dummy_num++;
					string dummyname = ss.str();
					newemap.insert(make_pair(dummyname,eweight));
					outDAG.nodeMap.insert(make_pair(dummyname,0));
					map<string, double> dummyedgemap ;
					dummyedgemap.insert(make_pair(dstr,eweight));
					outDAG.edgeMap.insert(make_pair(dummyname,dummyedgemap));
					
				}
				else{
					newemap.insert(make_pair(dstr,eweight));
				}
				
				
				
			}
			
			
		}
	  outDAG.edgeMap.insert(make_pair(node_id,newemap));
	}
	return outDAG;
}

vector<DAG> DAG::split(){
	vector<DAG> outVec;
	map<string, double> tmpnodeMap;
	map<string, double>::iterator nodeitr,neighboritr,addeditr;
	map<string, map<string, double> >::iterator edgeitr;
	tmpnodeMap = nodeMap;
	while(tmpnodeMap.size()>0){
		list<string > Q;		
		nodeitr = tmpnodeMap.begin();
		Q.push_back(nodeitr->first);
		

		map<string,double> addednodesMap;
		addednodesMap.insert(make_pair(nodeitr->first,0.0));
		while(Q.size()>0){
			
			string sstr = Q.front();
	 		Q.pop_front();
	 		
	 		//cout<< "add " <<sstr <<endl;
	 	 //sstr is the current node.
		 edgeitr = edgeMap.find(sstr);
		 map<string,double> neighborMap = edgeitr->second;
		 
		 for(neighboritr = neighborMap.begin();neighboritr !=neighborMap.end();neighboritr++){
		 		addeditr = addednodesMap.find(neighboritr->first);
		 		if(addeditr==addednodesMap.end()){
		 			
		 			Q.push_back(neighboritr->first);
		 			addednodesMap.insert(make_pair(neighboritr->first,0.0));
		 		}
		 }
		 
		}
		
		//pack DAG
		
		DAG newdag;
		for(addeditr = addednodesMap.begin();addeditr != addednodesMap.end();addeditr++){
			nodeitr = nodeMap.find(addeditr->first);
			newdag.nodeMap.insert(make_pair(nodeitr->first, nodeitr->second));
			tmpnodeMap.erase(nodeitr->first);
			edgeitr = edgeMap.find(addeditr->first);
			newdag.edgeMap.insert(make_pair(edgeitr->first, edgeitr->second));
		}
		
		//newdag.print();
		outVec.push_back(newdag);
		
		
	}
	return outVec;
}