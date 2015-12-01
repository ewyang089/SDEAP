#include "gene.h"
using namespace std;

SPGraph gene::build_from_GTF(std::string gtf_gname){

  SPGraph outGraph ;
	outGraph.junction =0;
	TOK tok;
	stringstream ss;
  TextIO tio;
  vector< vector<string> > gtfvec = tio.dlread(gtf_gname,"\t");
	gs.initGenomeSegment(gtf_gname,INT_MAX);

//add nodes and exon edges
	
	map<string, map<string, string> >::iterator exon_itr;
	for(exon_itr = gs.ExonMap.begin();exon_itr != gs.ExonMap.end() ; exon_itr++){
		
		if((exon_itr->second).size()>0){
			vector<string> tokvec = tok.Tokenize( exon_itr->first, "_" );
			ss.str("");
			
			ss<< tokvec[0] << "$"<<tokvec[1] ;
			outGraph.addNode(ss.str(), 0.0);
			
			
		}
	}
	
	
	
	
// add junction edge from GTF
	//parse all isoforms
	
	map<string, vector<string>  > IsoformMap;
	
	//init IsoformMap
	for(int i =0;i < gtfvec.size(); i++){
		map<string,string> gid_map = gs.parse_gid(gtfvec[i][8]);
		string t_id_str =  (gid_map.find("transcript_id"))->second ;
		map<string, vector<string>  >::iterator iso_itr = IsoformMap.find(t_id_str);
		if(gtfvec[i][2].compare("exon")==0){
			if( iso_itr == IsoformMap.end()){
				vector<string> tmpvec; 
				tmpvec.push_back(gtfvec[i][3]+"_"+gtfvec[i][4]);
			
				IsoformMap.insert(make_pair( t_id_str , tmpvec) );
			
			}
			else{
				(iso_itr->second).push_back(gtfvec[i][3]+"_"+gtfvec[i][4]);
			
			}
		}
	} 
	
	//find isoform ;
	
	map<string, vector<string>  >::iterator isomap_itr; 
	for(isomap_itr = IsoformMap.begin();isomap_itr != IsoformMap.end() ; isomap_itr++){
	
		vector<string> strvec = isomap_itr->second;
		map<string, vector<string> > exonseg_map;
		
		//add exon edges
		for(int i =0; i< strvec.size(); i++){
					vector<string> posstrvec = tok.Tokenize( strvec[i], "_");
					//cout << "exon : " << strvec[i] << endl;
					int pos_a = atoi(posstrvec[0].c_str());
					int pos_b = atoi(posstrvec[1].c_str());
					vector<string> nodeVec;
					map<string, string> addedMap;
					//search for segments to be linked
					for(int j = 0; j< gs.PosVec.size()-1 ; j++){
						if(pos_a <= gs.PosVec[j] && gs.PosVec[j] <= pos_b ){
							
							ss.str("");
							ss<<gs.PosVec[j] << "$"<< gs.PosVec[j];	
							if(outGraph.NodeMap.find(ss.str()) != outGraph.NodeMap.end()){
								if(addedMap.find(ss.str()) ==addedMap.end() ){
									nodeVec.push_back(ss.str());
									addedMap.insert(make_pair( ss.str(), ""));
								}
							}
							
							if(pos_a <= gs.PosVec[j+1] && gs.PosVec[j+1] <= pos_b ){
								ss.str("");
								ss<<gs.PosVec[j] << "$"<< gs.PosVec[j+1];	
								if(outGraph.NodeMap.find(ss.str()) != outGraph.NodeMap.end()){
									if(addedMap.find(ss.str()) ==addedMap.end() ){
										nodeVec.push_back(ss.str());
										addedMap.insert(make_pair( ss.str(), ""));
									}
								}
								ss.str("");
								ss<<gs.PosVec[j+1] << "$"<< gs.PosVec[j+1];	
								if(outGraph.NodeMap.find(ss.str()) != outGraph.NodeMap.end()){
									if(addedMap.find(ss.str()) ==addedMap.end() ){
										nodeVec.push_back(ss.str());
										addedMap.insert(make_pair( ss.str(), ""));
									}
								}
								
							}
							
						}
						
						
					}			
					
					//linked segments
					for(int j =0 ; j< nodeVec.size()-1; j++){
						//cout << nodeVec[j] << "->"<< nodeVec[j+1] <<endl;
						ss.str("");
						ss << nodeVec[j]<< "_"<< nodeVec[j+1];
						if(outGraph.EdgeMap.find(ss.str()) == outGraph.EdgeMap.end()){
							outGraph.addEdge(nodeVec[j], nodeVec[j+1],  "E",0 ,0 );
						}
					}
					exonseg_map.insert(make_pair(strvec[i] , nodeVec));
					
		} 
		
		//add junction edges
		for( int i=0; i< strvec.size()-1; i++){
			
			vector<string> nodevec1 = exonseg_map.find(strvec[i])->second;
			vector<string> nodevec2 = exonseg_map.find(strvec[i+1])->second;
			ss.str("");
			string n1_str = nodevec1[nodevec1.size()-1];
			string n2_str = nodevec2[0];
			
			ss<< n1_str<< "_"<< n2_str;
			
			if(outGraph.EdgeMap.find(ss.str()) == outGraph.EdgeMap.end()){
				outGraph.addEdge(n1_str, n2_str, "J",0 ,0 );	
			}
			
			
		}
		
		//cout << endl;
	}	 
	
	
	//addd s and t 
	
	
	map<string, Node >::iterator n_itr; 
	vector<string> in_vec,out_vec ;	



	for(isomap_itr = IsoformMap.begin();isomap_itr != IsoformMap.end() ; isomap_itr++){
		  vector<string> exonvec=  isomap_itr -> second;
			//find the head node of the isoform
			vector<string> tmpvec = tok.Tokenize(exonvec[0], "_");
			
			string idx_str = gs.find_interval(atoi(tmpvec[0].c_str()));
			ss.str("");	
			tmpvec = tok.Tokenize(idx_str, "_");
			ss<< gs.PosVec[atoi(tmpvec[0].c_str())]<< "$"<< gs.PosVec[atoi(tmpvec[1].c_str())];
			in_vec.push_back(ss.str());
			//find the end node of the isoform 
			
			tmpvec = tok.Tokenize(exonvec[exonvec.size()-1], "_");
			idx_str = gs.find_interval(atoi(tmpvec[1].c_str()));
			ss.str("");	
			tmpvec = tok.Tokenize(idx_str, "_");
			ss<< gs.PosVec[atoi(tmpvec[0].c_str())]<< "$"<< gs.PosVec[atoi(tmpvec[1].c_str())];
			out_vec.push_back(ss.str());
			
				
	}



outGraph.addNode("S", 0.0);
	outGraph.addNode("T", 0.0);

	for(int i =0;i< in_vec.size() ; i++){
		
		outGraph.addEdge("S", in_vec[i], "J",0 ,0 );
	}
	
	for(int i =0;i< out_vec.size() ; i++){
		outGraph.addEdge(out_vec[i],"T", "J",0 ,0 );
	}
	
	return outGraph;
	
}


void gene::map_one_read_bed(std::string read1_str, double scale){
	string str =read1_str;
	TOK tok ;
	stringstream ss;
	vector<string> tokvec = tok.Tokenize(str, "\t");
			
	int read_spos = atoi(tokvec[1].c_str() );
	read_spos= read_spos+1; // BED is 0 based
	vector<string> sizevec = tok.Tokenize(tokvec[10], ",");
	vector<string> offsetvec = tok.Tokenize(tokvec[11], ",");
			
	int len =0;
	for(int i =0; i < sizevec.size() ; i++){
		len += atoi(sizevec[i].c_str()); 
				
	}
			
	for(int i =0;i<offsetvec.size();i++){
				
				int size = atoi(sizevec[i].c_str());
				int spos = read_spos + atoi(offsetvec[i].c_str());
				int epos = spos + size -1 ;
				//cout << "seg pos :"<< spos<< ","<< epos << endl;
				//find segs 
				vector<string> intvec = gs.get_intervals_by_pos(spos,epos);
				for(int j =0; j< intvec.size();j++){
					vector<string> tmpvec = tok.Tokenize(intvec[j], "_");
					int int_posa= atoi(tmpvec[0].c_str()); 
					int int_posb= atoi(tmpvec[1].c_str());
					//cout <<"chk int: "<< intvec[j] << endl;
					if( int_posa <= spos && epos <= int_posb){
					// this read is contained in the interval. 	
						ss.str("");
						ss<< int_posa<< "$"<< int_posb;
						map<string,Node>::iterator node_itr =  spg.NodeMap.find(ss.str());
						if(node_itr != spg.NodeMap.end()){
							(node_itr->second).weight = (node_itr->second).weight+ ((size/(double)len)* scale);
						}
					}
					else{
						if(int_posa < spos ){
							
							double ov =  int_posb - spos +1;
							ov = ov/(double)len;
							ss.str("");
							ss<< int_posa<< "$"<< int_posb;
							//cout << "kuku:" << ss.str()<< ","<< ov<< endl;
							map<string,Node>::iterator node_itr =  spg.NodeMap.find(ss.str());
							if(node_itr != spg.NodeMap.end()){
								(node_itr->second).weight = (node_itr->second).weight+ (ov* scale);
							}	
							
						}else if(int_posb > epos){
							
							double ov = epos - int_posa   +1;
							ov = ov/(double)len;
							ss.str("");
							ss<< int_posa<< "$"<< int_posb;
							map<string,Node>::iterator node_itr =  spg.NodeMap.find(ss.str());
							if(node_itr != spg.NodeMap.end()){
								(node_itr->second).weight = (node_itr->second).weight+ (ov* scale);
							}	
						}
						else{
							double ov = int_posb - int_posa   +1;
							ov = ov/(double)len;
							ss.str("");
							ss<< int_posa<< "$"<< int_posb;
							map<string,Node>::iterator node_itr =  spg.NodeMap.find(ss.str());
							if(node_itr != spg.NodeMap.end()){
								(node_itr->second).weight = (node_itr->second).weight+ (ov* scale);
							}
							
						}
						
					}
					
				}
				 
				
	}
	
}

void gene::create_coverage(std::string gtf_gname, std::string bed_fname){
	gs.initGenomeSegment(gtf_gname,INT_MAX);
}

void gene::map_one_read_bed_toSPG(std::string read1_str, double scale){
	string str =read1_str;
	TOK tok ;
	stringstream ss;
	vector<string> tokvec = tok.Tokenize(str, "\t");
			
	int read_spos = atoi(tokvec[1].c_str() );
	read_spos= read_spos+1; // BED is 0 based
	vector<string> sizevec = tok.Tokenize(tokvec[10], ",");
	vector<string> offsetvec = tok.Tokenize(tokvec[11], ",");
			
	int len =0;
	for(int i =0; i < sizevec.size() ; i++){
		len += atoi(sizevec[i].c_str()); 
				
	}
		
	if(offsetvec.size()<=2){
		
		if(offsetvec.size() ==1)
		{
			
			int size = atoi(sizevec[0].c_str());
			int spos = read_spos + atoi(offsetvec[0].c_str());
			int epos = spos + size -1 ;
			
			vector<string> intvec = gs.get_intervals_by_pos(spos,epos);
			
			if(intvec.size()==1){
				//fully in a seg
				vector<string> tmpvec = tok.Tokenize(intvec[0], "_");
					int int_posa= atoi(tmpvec[0].c_str()); 
					int int_posb= atoi(tmpvec[1].c_str());
										if( int_posa <= spos && epos <= int_posb){
					// this read is contained in the interval. 	
						ss.str("");
						ss<< int_posa<< "$"<< int_posb;
						map<string,Node>::iterator node_itr =  spg.NodeMap.find(ss.str());
						if(node_itr != spg.NodeMap.end()){
							
							(node_itr->second).weight = (node_itr->second).weight+ (1.0* scale);
							
						}
					}
			}
			
			if(intvec.size()==2){
				//exon link
				string nodeA;
				string nodeB;
				
				for(int j =0; j <2;j++){
				vector<string> tmpvec = tok.Tokenize(intvec[j], "_");
					int int_posa= atoi(tmpvec[0].c_str()); 
					int int_posb= atoi(tmpvec[1].c_str());
					ss.str("");
					ss<< int_posa<< "$"<< int_posb;
					if(j==0)
					{
						nodeA = ss.str(); 
					}
					if(j==1)
					{
						nodeB = ss.str(); 
					}
				}
				
				map<string,Edge>::iterator edge_itr =  spg.EdgeMap.find( nodeA+"_"+nodeB );
			  if(edge_itr != spg.EdgeMap.end()){
			  			
							(edge_itr->second).weight = (edge_itr->second).weight+ (1.0* scale);
							
				}
				
			}	
		}
		if(offsetvec.size() ==2){
			//junction
			
			
			int size = atoi(sizevec[0].c_str())+atoi(sizevec[1].c_str());
			int spos1 = read_spos + atoi(offsetvec[0].c_str());
			int epos1 = spos1 + atoi(sizevec[0].c_str()) -1 ;
			int spos2 = read_spos + atoi(offsetvec[1].c_str());
			int epos2 = spos2 + atoi(sizevec[1].c_str()) -1 ;
			string JnodeA,JnodeB;
			JnodeA = JnodeB ="";
			//part A
			vector<string> intvec = gs.get_intervals_by_pos(spos1,epos1);
			if(intvec.size() >0 ){
				vector<string> tmpvec = tok.Tokenize(intvec[intvec.size()-1], "_");
				int int_posa= atoi(tmpvec[0].c_str()); 
				int int_posb= atoi(tmpvec[1].c_str());
						// this read is contained in the interval. 	
				ss.str("");
				ss<< int_posa<< "$"<< int_posb;
				JnodeA= ss.str();
			}	
			//part B
			intvec = gs.get_intervals_by_pos(spos2,epos2);
			
			if(intvec.size()>0){
				vector<string> tmpvec = tok.Tokenize(intvec[0], "_");
				int int_posa= atoi(tmpvec[0].c_str()); 
				int int_posb= atoi(tmpvec[1].c_str());
				ss.str("");
				ss<< int_posa<< "$"<< int_posb;
				JnodeB= ss.str();
			}	
			//find edge
			
				map<string,Edge>::iterator edge_itr =  spg.EdgeMap.find( JnodeA+"_"+JnodeB );
			  if(edge_itr != spg.EdgeMap.end()){
			  	
							(edge_itr->second).weight = (edge_itr->second).weight+ (1.0* scale);
							
				}
			
		}
	}		
	
		
}



void gene::map_one_frag_bed(std::string read1_str , std::string read2_str){
		 		map_one_read_bed(read1_str, 0.5 );
		 		map_one_read_bed(read2_str, 0.5 );
	
}



void gene::map_one_bed_file(std::string fname,double scale){
	stringstream ss;
	fstream filestr;
	
	
	filestr.open(fname.c_str(), ios_base::in);
	
	string str;	 
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
	 		
	 		map_one_read_bed(str, scale);
	 	
	 	}
	}

	filestr.close();
}


void gene::map_one_bed_file_toSPG(std::string fname, double scale){
	stringstream ss;
	fstream filestr;
	
	
	filestr.open(fname.c_str(), ios_base::in);
		
	string str;	 
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
	 		//cout << "map:" << str <<endl;
	 		map_one_read_bed_toSPG(str, scale);
	 		//cout << "map:" << str <<endl;
	 	}
	}
	
	//spg._debug();
	filestr.close();
}

int  gene::multiGene(std::string gtf_gname){
	int outD=0;
	spg = build_from_GTF(gtf_gname);
	
	D.decomposeG(spg);
	//postoder
	outD = D.ASM_Map.size();
	
	return outD;
}

inline bool exists_test (const std::string& name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}

std::vector<std::string> gene::init(std::string gtf_gname,std::string bed_fname,std::string gname,double scale){
	vector<string> outVec;
	if(exists_test(bed_fname.c_str()) && exists_test(gtf_gname.c_str())){
	spg = build_from_GTF(gtf_gname);
	spg.gname = gname;
	map_one_bed_file(bed_fname,scale);
	
	stringstream ss;
	
	ss.precision(15);
	for(map<string, Node>::iterator node_itr = spg.NodeMap.begin() ; node_itr != spg.NodeMap.end();node_itr++){
		ss.str("");
		Node & node_u = node_itr->second;
		double rc = node_u.weight ;
		ss<< node_itr->first << "\t" << rc; 
		outVec.push_back(ss.str()) ;
	}
	
	}
	return outVec;
}



std::vector<std::string> gene::init_SPG(std::string gtf_gname,std::string bed_fname,std::string gname,string tmpdir, int reads_len, double scale){
	
	TextIO tio;
	TOK tok;
	stringstream ss;
	
	vector<string> weightvec,gene_fpkm_vec ;
	
	if(exists_test(bed_fname.c_str()) && exists_test(gtf_gname.c_str())){
	spg = build_from_GTF(gtf_gname);
	
	//init add pseudo count to 
/*	for(map<string, Edge>::iterator edge_itr = spg.EdgeMap.begin() ; edge_itr != spg.EdgeMap.end() ; edge_itr++){
		Edge & edge = edge_itr->second;
		edge.weight = scale;
	}
*/	
	//
  
	spg.gname = gname;
	map_one_bed_file_toSPG(bed_fname,scale);
	//map_one_bed_file(bed_fname,scale);
	spg._debug();
	D.decomposeG(spg);
	

	int i =0;
	//evaluate every ASM in the gene
	for(map<std::string, SPGraph>::iterator asm_itr =  D.ASM_Map.begin(); asm_itr !=  D.ASM_Map.end() ;asm_itr ++){
		
		cout << asm_itr->first << "=============================================================" << endl;
		ss.str("");
		(asm_itr->second)._debug();
		
		//if((asm_itr->first).compare("S~T") != 0){
			vector<vector<string> > mtxvec =  D.toMTX_SPG(asm_itr->second,reads_len);
			if(mtxvec.size()>0 ){
				ss<< gname<<"."<<++i<<".mtx";
				
				tio.dlwrite( tmpdir+"/"+ss.str(), mtxvec,"\t");
			
				string cmd = "Rscript OPT.r "+ tmpdir+"/"+ss.str();
				system(cmd.c_str());
				cmd = "rm -rf "+ tmpdir+"/"+ss.str();
				//system(cmd.c_str());
			}
			weightvec.push_back(tmpdir+"/"+ss.str()+".weight");
		//}
		/*else{
			
			  vector<string> outVec;
			  double len =0;
			  double total_rc =0;
			  ss.str("");
				for(map<string, Node>::iterator node_itr = spg.NodeMap.begin() ; node_itr != spg.NodeMap.end();node_itr++){
					
					Node & node_u = node_itr->second;
					double rc = node_u.weight ;
					
					vector<string> namevec =  tok.Tokenize(node_itr->first,"$");
					int posa =atof(namevec[0].c_str()); 
					int posb =atof(namevec[1].c_str());
					double exon_len = posb - posa +1;
					ss<< "N:"<<node_itr->first <<","; 
					len+=exon_len;
					total_rc+= rc;
				}
				
				for(map<string, Edge>::iterator edge_itr = spg.EdgeMap.begin() ; edge_itr != spg.EdgeMap.end();edge_itr++){
					
					Edge & edge_u = edge_itr->second;
					double rc = edge_u.weight ;
					total_rc+= rc;
				}
				
				ss<< "\t" << total_rc*1000 / (double) len;
				outVec.push_back(ss.str()) ;
				ss.str("");
				ss<< gname<<"."<<++i<<".mtx";
				
				tio.dlwrite(tmpdir+"/"+ss.str()+".weight", outVec);
				//cout << "writing :" << tmpdir+"/"+ss.str()+".weight"<< endl;
				weightvec.push_back(tmpdir+"/"+ss.str()+".weight");
		}*/
		
		
	}
	//compile rpkms of all paths
  
	
	for(int i =0; i < weightvec.size() ; i++){
		//cout << weightvec[i]<< endl;
		vector<vector<string> > asmvec =  tio.dlread(weightvec[i],"*");
		
		for(int j=0 ; j<asmvec.size(); j++){
			gene_fpkm_vec.push_back(asmvec[j][0]);
			
		}
		string cmd = "rm -rf "+ weightvec[i];
		//system(cmd.c_str());
	}
	
	//tio.dlwrite(tmpdir+"/"+gname+".fpkm", gene_fpkm_vec);
	}	
	
	return gene_fpkm_vec;
	
	
}

std::vector<std::string> gene::init_Hybrid(std::string gtf_gname,std::string bed_fname,std::string gname,string tmpdir, int reads_len, double scale){
	
	TextIO tio;
	TOK tok;
	stringstream ss;
	
	vector<string> weightvec,gene_fpkm_vec ;
			
	if(exists_test(bed_fname.c_str()) && exists_test(gtf_gname.c_str())){
	spg = build_from_GTF(gtf_gname);

 	
	spg.gname = gname;
	map_one_bed_file_toSPG(bed_fname,scale);

	D.decomposeG(spg);
	
	
	int ASM_n =0;
	//evaluate every ASM in the gene
	for(map<std::string, SPGraph>::iterator asm_itr =  D.ASM_Map.begin(); asm_itr !=  D.ASM_Map.end() ;asm_itr ++){
		
		
		SPGraph & ASM_G = (asm_itr->second);
		//cout << asm_itr->first << "=============================================================" << endl;
		//(asm_itr->second)._debug();
		
		vector<string> pathvec = ASM_G.list_Paths();
		//cout << ASM_n+1<< ":" << pathvec.size()<< endl;
		if(pathvec.size()<=3){
		
			vector<vector<string> > mtxvec =  D.toMTX_SPG(ASM_G,reads_len);
		
			if(mtxvec.size()>0 ){
				ss.str("");
				ss<< gname<<"."<<++ASM_n<<".mtx";
				
				tio.dlwrite( tmpdir+"/"+ss.str(), mtxvec,"\t");
				weightvec.push_back(ss.str()+".weight");
				
				string cmd = "Rscript OPT.r "+ tmpdir+"/"+ss.str();
				system(cmd.c_str());
				cmd = "rm -rf "+ tmpdir+"/"+ss.str();
				system(cmd.c_str());
				
				
			}
			
		}
		else{
			ASM_n++;
		}
				
	}
	//compile rpkms of all paths
  

	}	
	
	return weightvec;
	
	
}

std::vector<std::string> gene::init_Hybrid_RC(std::string gtf_gname,std::string bed_fname,std::string gname,string tmpdir, int reads_len, double scale){
	
	TextIO tio;
	TOK tok;
	stringstream ss;
	
	vector<string> weightvec,gene_fpkm_vec ;
	
	if(exists_test(bed_fname.c_str()) && exists_test(gtf_gname.c_str())){
	spg = build_from_GTF(gtf_gname);
	//cout << bed_fname << endl;
	
	//init add pseudo count to 
	
	//
  
	spg.gname = gname;
	map_one_bed_file(bed_fname,scale);
	/*
	double sum = 0;
	for(map<string, Node>::iterator node_itr = spg.NodeMap.begin() ; node_itr != spg.NodeMap.end() ; node_itr++){
	
		Node & node = node_itr->second;
		sum += node.weight;
	}
	
	for(map<string, Edge>::iterator edge_itr = spg.EdgeMap.begin() ; edge_itr != spg.EdgeMap.end() ; edge_itr++){
		Edge & edge = edge_itr->second;
		sum+= edge.weight ;
	}
	*/
	//cout << sum/(double)scale << endl;
	
	//map_one_bed_file(bed_fname,scale);
	D.decomposeG(spg);
	
	
	int ASM_n =0;
	//evaluate every ASM in the gene
	for(map<std::string, SPGraph>::iterator asm_itr =  D.ASM_Map.begin(); asm_itr !=  D.ASM_Map.end() ;asm_itr ++){
		//cout << asm_itr->first << "============================" << endl;
		ss.str("");
		SPGraph & ASM_G = (asm_itr->second);
	
		vector<string> pathvec = ASM_G.list_Paths();
//		cout << ASM_n+1<< ":" << pathvec.size()<< endl;
		if(pathvec.size()<=3){
		
			ASM_n++;
			
		}
		else{
				vector<string> ASM_RCVec;
			//break it down into (collasped) exons
				map<string, string> KeyNodeMap ;
				map<string,string> addedmap;
				for(map<string, Node>::iterator n_itr = ASM_G.NodeMap.begin(); n_itr != ASM_G.NodeMap.end() ; n_itr++){
					Node & node_u = n_itr->second;
					if(node_u.in_d> 1 || node_u.out_d> 1){
						KeyNodeMap.insert( make_pair(n_itr->first, ""));
					}
					
				}
				
				for(map<string, Node>::iterator n_itr = ASM_G.NodeMap.begin(); n_itr != ASM_G.NodeMap.end() ; n_itr++){
					Node & node_u = n_itr->second;
					if((n_itr->first).compare("S")!=0 && (n_itr->first).compare("T")!=0){
					if(KeyNodeMap.find(n_itr->first) != KeyNodeMap.end()){
							ss.str("");
							vector<string> posvec = tok.Tokenize(n_itr->first,"$");
							int len = atof(posvec[1].c_str()) - atof(posvec[0].c_str()) +1;
							ss<< n_itr->first << ",\t" << len<< "\t"<< node_u.weight<< "\t" << 1000*node_u.weight/len;
							string feature_name = n_itr->first+",";
							if(addedmap.find(feature_name) == addedmap.end()){
								ASM_RCVec.push_back(ss.str());
								addedmap.insert(make_pair(feature_name , ""));
								
							} 
						}
					}
					
				}
				
				vector<string> ExonVec = ASM_G.list_UniPaths();
				
					//collaspe exons
				for(int i =0 ; i < ExonVec.size() ; i++){
				
						
						vector<string> nodevec = tok.Tokenize(ExonVec[i],",");
						if(nodevec.size()>1){
						ss.str("");
						double rc =0;
						double len =0;
						string feature_name ="";
						for(int j=0; j<nodevec.size(); j++){
							
							
							if(nodevec[j].compare("S") != 0 && nodevec[j].compare("T") != 0){
								if(KeyNodeMap.find(nodevec[j]) == KeyNodeMap.end()){
									ss<<nodevec[j] <<",";
									Node & node_u = ASM_G.NodeMap.find(nodevec[j])->second;
								
									rc +=  node_u.weight;
									vector<string> posvec = tok.Tokenize(nodevec[j],"$");
									len+= atof(posvec[1].c_str()) - atof(posvec[0].c_str()) +1;
								}
							
							}
						}
						feature_name = ss.str();
						if(feature_name.length() >0){
							ss<<"\t" <<len << "\t" << rc << "\t"<< 1000*rc/len;
							
							if(addedmap.find(feature_name) == addedmap.end()){
								ASM_RCVec.push_back(ss.str());
								addedmap.insert(make_pair(feature_name , ""));
							}
						}
						
						
						}
					
				}
				ss.str("");
				ss<< gname<<"."<<++ASM_n;
				tio.dlwrite( tmpdir+"/"+ss.str()+".rc", ASM_RCVec);
				weightvec.push_back(ss.str()+".rc");
		}
		
		
	}
	//compile rpkms of all paths
  
	/*
	for(int i =0; i < weightvec.size() ; i++){
		//cout << weightvec[i]<< endl;
		vector<vector<string> > asmvec =  tio.dlread(weightvec[i],"*");
		
		for(int j=0 ; j<asmvec.size(); j++){
			gene_fpkm_vec.push_back(asmvec[j][0]);
			
		}
		string cmd = "rm -rf "+ weightvec[i];
		//system(cmd.c_str());
	}
	*/
	//tio.dlwrite(tmpdir+"/"+gname+".fpkm", gene_fpkm_vec);
	}	
	
	return weightvec;
	
	
}



void gene::init_SPG_exon(std::string gtf_gname,std::string bed_fname,std::string gname,string tmpdir, int reads_len, double scale){
	
	TextIO tio;
	TOK tok;
	stringstream ss;
	
	vector<string> weightvec,gene_fpkm_vec ;
	
	if(exists_test(bed_fname.c_str()) && exists_test(gtf_gname.c_str())){
		spg = build_from_GTF(gtf_gname);
  	//spg._debug();
		spg.gname = gname;
	//map_one_bed_file_toSPG(bed_fname,scale);
  	map_one_bed_file(bed_fname,scale);
		D.decomposeG(spg);
	
		
	//evaluate every ASM in the gene
		
		vector<string> ASM_RCVec, ST_RCVec; 
		map<string,string> addedmap;
		for(map<std::string, SPGraph>::iterator asm_itr =  D.ASM_Map.begin(); asm_itr !=  D.ASM_Map.end() ;asm_itr ++){
			string asm_name = asm_itr->first;
			//cout << asm_name << "=============================================================" << endl;
			SPGraph & ASM_G = asm_itr->second;
			ss.str("");
			
			//ASM_G._debug();
			
				map<string, string> KeyNodeMap ;
				for(map<string, Node>::iterator n_itr = ASM_G.NodeMap.begin(); n_itr != ASM_G.NodeMap.end() ; n_itr++){
					Node & node_u = n_itr->second;
					if(node_u.in_d> 1 || node_u.out_d> 1){
						KeyNodeMap.insert( make_pair(n_itr->first, ""));
					}
					
				}
				
				for(map<string, Node>::iterator n_itr = ASM_G.NodeMap.begin(); n_itr != ASM_G.NodeMap.end() ; n_itr++){
					Node & node_u = n_itr->second;
					if((n_itr->first).compare("S")!=0 && (n_itr->first).compare("T")!=0){
					if(KeyNodeMap.find(n_itr->first) != KeyNodeMap.end()){
							ss.str("");
							vector<string> posvec = tok.Tokenize(n_itr->first,"$");
							int len = atof(posvec[1].c_str()) - atof(posvec[0].c_str()) +1;
							ss<< n_itr->first << ",\t" << len<< "\t"<< node_u.weight<< "\t" << 1000*node_u.weight/len;
							string feature_name = n_itr->first+",";
							if(addedmap.find(feature_name) == addedmap.end()){
								if(asm_name.compare("S~T") !=0){
										ASM_RCVec.push_back(ss.str());
								}else{
										ST_RCVec.push_back(ss.str());
								}
								addedmap.insert(make_pair(feature_name , ""));
								
							} 
						}
					}
					
				}
				
				vector<string> ExonVec = ASM_G.list_UniPaths();
				
					//collaspe exons
				for(int i =0 ; i < ExonVec.size() ; i++){
				
						
						vector<string> nodevec = tok.Tokenize(ExonVec[i],",");
						if(nodevec.size()>1){
						ss.str("");
						double rc =0;
						double len =0;
						string feature_name ="";
						for(int j=0; j<nodevec.size(); j++){
							
							
							if(nodevec[j].compare("S") != 0 && nodevec[j].compare("T") != 0){
								if(KeyNodeMap.find(nodevec[j]) == KeyNodeMap.end()){
									ss<<nodevec[j] <<",";
									Node & node_u = ASM_G.NodeMap.find(nodevec[j])->second;
								
									rc +=  node_u.weight;
									vector<string> posvec = tok.Tokenize(nodevec[j],"$");
									len+= atof(posvec[1].c_str()) - atof(posvec[0].c_str()) +1;
								}
							
							}
						}
						feature_name = ss.str();
						if(feature_name.length() >0){
							ss<<"\t" <<len << "\t" << rc << "\t"<< 1000*rc/len;
							
							if(addedmap.find(feature_name) == addedmap.end()){
								if(asm_name.compare("S~T") !=0){
									ASM_RCVec.push_back(ss.str());
								}else{
									ST_RCVec.push_back(ss.str());
								}
								addedmap.insert(make_pair(feature_name , ""));
							}
						}
						
						
						}
					
				}
		
		}
		
	 tio.dlwrite(tmpdir+"/"+gname+".ASM_RC", ASM_RCVec);
	 tio.dlwrite(tmpdir+"/"+gname+".ST_RC", ST_RCVec);
	}	
	

	
}