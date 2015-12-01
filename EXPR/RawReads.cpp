#include "RawReads.h"



using namespace std;


void RawReads::sep_gtfs(std::string fnamestr, std::string chrname){
	TextIO tio;
	TOK tok;
	map<string, vector<string> > ChrmReadMap;
	stringstream ss;
	

	
	fstream filestr;

	filestr.open(fnamestr.c_str(), ios_base::in);
		
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			vector<string> tmpvec = tok.Tokenize(str,"\t" );
			if(tmpvec[0].compare(chrname)==0){
				cout << str << endl;
			}
		}
	}
	
	filestr.close();
	
	
	
}


void RawReads::sep_gtfs_by_chrs(std::string fnamestr,std::string out_dir){
	TextIO tio;
	TOK tok;
	map< string, vector<string> > gtfmap;
	stringstream ss;

	fstream filestr;

	filestr.open(fnamestr.c_str(), ios_base::in);
	
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			vector<string> tmpvec = tok.Tokenize(str,"\t" );
			if( gtfmap.find(tmpvec[0]) == gtfmap.end()){
				vector<string> tokvec;
				//cout << tmpvec[1] << endl;
				tokvec.push_back(str);
				gtfmap.insert(make_pair(tmpvec[0],tokvec));
			}
			else{
				
				(gtfmap.find(tmpvec[0])->second).push_back(str);
				
			}
		}
	}
	
	filestr.close();
	
	for(map< string, vector<string> >::iterator mapitr = gtfmap.begin();mapitr != gtfmap.end() ; mapitr++){
		tio.dlwrite(out_dir+"/"+mapitr->first,  mapitr->second );
	}
	
}

void RawReads::sep_reads_by_chrs(string out_dir , std::string fnamestr){
		TextIO tio;
	TOK tok;
	map< string, vector<string> > gtfmap;
	stringstream ss;

	fstream filestr;

	filestr.open(fnamestr.c_str(), ios_base::in);
	
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			vector<string> tmpvec = tok.Tokenize(str,"\t" );
			if( gtfmap.find(tmpvec[0]) == gtfmap.end()){
				vector<string> tokvec;
				//cout << tmpvec[1] << endl;
				tokvec.push_back(str);
				gtfmap.insert(make_pair(tmpvec[0],tokvec));
			}
			else{
				
				(gtfmap.find(tmpvec[0])->second).push_back(str);
				
			}
		}
	}
	
	filestr.close();
	
	for(map< string, vector<string> >::iterator mapitr = gtfmap.begin();mapitr != gtfmap.end() ; mapitr++){
		tio.dlwrite(out_dir+"/"+mapitr->first,  mapitr->second );
	}
	
	
}

void RawReads::sep_reads_by_genes(std::string out_dir , std::string gtfname, std::string bedname){
	GenomeSegment gs;
	gs.initGenomeSegment(gtfname, INT_MAX, "gene_id");

	TOK tok;
	fstream filestr;
	map< string, vector<string> > RMap;
	TextIO tio;
	filestr.open(bedname.c_str(), ios_base::in);
		
	string str;	 
	while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 	
	 		vector<string> tokvec = tok.Tokenize( str,"\t");
	 		int spos = atoi(tokvec[1].c_str());
	 		spos = spos+1;
	 		int epos = atoi(tokvec[2].c_str());	
	 		vector<string> segvec =  gs.get_intervals_by_pos(spos,epos);
	 		map<string, string> AddedMap;
	 		for(int i =0; i<segvec.size(); i++){
	 			map<string, map<string,string> >::iterator exon_itr = gs.ExonMap.find(segvec[i]);
	 				if(exon_itr!=gs.ExonMap.end()){
	 					map<string,string> &gmap = exon_itr->second;
	 					
	 					for(map<string,string>::iterator g_itr = gmap.begin(); g_itr != gmap.end(); g_itr++){
	 						if(AddedMap.find(g_itr->first) == AddedMap.end()){
	 							map< string , vector<string> >::iterator R_itr = RMap.find(g_itr->first);
	 							if( R_itr!= RMap.end() ){
	 								(R_itr->second).push_back(str);
	 								AddedMap.insert(make_pair(g_itr->first,""));
	 							}
	 							else{
	 								vector<string> tmpvec;
	 								tmpvec.push_back(str);
	 								RMap.insert(make_pair( g_itr->first, tmpvec));
	 								AddedMap.insert(make_pair(g_itr->first,""));
	 							}
	 						}
	 					}
	 				}
	 		}
		}
	}
	filestr.close();
	int all_rc = 0 ;
	for(map<string, vector<string> >::iterator r_itr = RMap.begin(); r_itr != RMap.end(); r_itr++){
		
		tio.dlwrite( out_dir+"/"+ r_itr->first,r_itr->second );
		all_rc += (r_itr->second).size();
	}
	
	cout << all_rc << endl;
}

void RawReads::paired_reads_by_chrs(std::string out_dir ){
	stringstream ss;
	TextIO tio;
	TOK tok;
	for(int i=1;i<=20;i++){
		map<string , vector<string> > PairedMap; 
		
		ss.str("");
		if(i < 20){
			ss<< i;
		}if(i==20){
			ss<<"X" ;
		}
		
		string chrmname = "chr" +ss.str();
		cout <<out_dir+"/"+chrmname << endl;
		
		vector<string> bedvec = tio.dlread_vector(out_dir+"/"+chrmname);
		
		for(int j =0 ; j< bedvec.size(); j++){
			vector<string> tmpvec = tok.Tokenize(bedvec[j],"\t");
			string bedname = tmpvec[3];
			tmpvec = tok.Tokenize(bedname,"/");
			string readname = tmpvec[0];
			map<string , vector<string> >::iterator pitr = PairedMap.find(readname);
			if(pitr == PairedMap.end()){
				vector<string> strvec ;
				strvec.push_back(bedvec[j]);
				PairedMap.insert( make_pair( readname, strvec ) );
				
			}
			else{
				(pitr->second).push_back(bedvec[j]);
			}
		}
		vector<string> outvec;
  	for(  map<string , vector<string> >::iterator mapitr = PairedMap.begin(); mapitr != PairedMap.end(); mapitr ++ ){
  		vector<string> pairs =(mapitr->second); 
  		if(pairs.size() ==2 ){
  			outvec.push_back(pairs[0]);
  			outvec.push_back(pairs[1]);
  		}
  	}
  	
  	tio.dlwrite(out_dir+"/"+chrmname+"_paired.bed", outvec);
  	
	}
}
	
	
void RawReads::sep_reads(std::string fnamestr, std::string chrname){
	TextIO tio;
	TOK tok;
	map<string, vector<string> > ChrmReadMap;
	stringstream ss;
	

	
	fstream filestr;

	filestr.open(fnamestr.c_str(), ios_base::in);
		
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			vector<string> tmpvec = tok.Tokenize(str,"\t" );
			if(tmpvec[0].compare(chrname)==0){
				cout << str << endl;
			}
		}
	}
	
	filestr.close();
	
	
	
}
std::string RawReads::reverse_string(std::string in_str){
	stringstream ss; 
	ss.str("");
	for(int i =0; i< in_str.length(); i++){
		char cc = in_str[in_str.length()-1-i];
		if(cc == 'A')
		{
			ss<<'T';
		}
		if(cc == 'T')
		{
			ss<<'A';
		}	
		if(cc == 'C')
		{
			ss<<'G';
		}
		if(cc == 'G')
		{
			ss<<'C';
		}
	}
	return ss.str();
}

string readFasta(string fname){
	string outstr=""; 	
	fstream filestr;
	
	string str;	 
	filestr.open(fname.c_str(), ios_base::in);
	getline(filestr,str);
	
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			
			outstr+=str;
			
		}
	}
	
	filestr.close();
	
	
	return outstr;
}

void RawReads::simulated_bed2Sam(std::string infname, std::string outfname){
	
//	cout << chrmap.find(chr_dir[0])->second<< endl;
	TextIO tio;
	TOK tok;
	stringstream ss;
	fstream filestr;
	vector<string> SamVec;
	
	filestr.open(infname.c_str(), ios_base::in);
		
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			ss.str("");
			vector<string> bedvec = tok.Tokenize(str , "\t");
			vector<string> sizevec = tok.Tokenize(bedvec[10] , ",");
			vector<string> offsetvec = tok.Tokenize(bedvec[11] , ",");
			int p_sum =0;
			for(int i =0 ; i<sizevec.size(); i++){
				//M
				int size = atoi(sizevec[i].c_str());
			  p_sum += size;
			  ss<< size	<<"M";
				if(i<sizevec.size()-1){
					int offset = atoi(offsetvec[i+1].c_str());
					ss<< offset-p_sum	<<"N";
					p_sum += offset-p_sum ; 
				}
			}
			string cigar = ss.str();
			ss.str("");
			ss<< bedvec[3] <<"\t" << 0<<"\t" << bedvec[0] << "\t"<< atoi(bedvec[1].c_str())+1 << "\t255\t"	<<cigar<<"\t*\t0\t0\t*\t*\t" << "NM:i:0 NH:i:1";
			if(sizevec.size()>1){
				ss<<" XS:A:+";
			}  
			SamVec.push_back(ss.str()); 
		}
	}
	
	tio.dlwrite( outfname,SamVec );
	
}


	
void RawReads::bed2Fasta_paired(string fname ,vector<string> chr_dir){

	map<string, string> readmap;		
	vector<string> PosReadsVec;

	vector<string> outvec ; 
	
	map<string,string> chrmap;
	//read chr strings
	for(int i =0 ; i< chr_dir.size();i++)
	{
		cout << "/rhome/ywyang/bigdata/RefSeq/hg38/chroms/" +chr_dir[i]+".fa" <<endl;
		string myString = readFasta("/rhome/ywyang/bigdata/RefSeq/hg38/chroms/" +chr_dir[i]+".fa");
	  
		std::transform(myString.begin(), myString.end(), myString.begin(), ::toupper);
    chrmap.insert(make_pair(chr_dir[i], myString));
    
	}
	
//	cout << chrmap.find(chr_dir[0])->second<< endl;
	TOK tok;
	stringstream ss;
	fstream filestr;
	
	
	filestr.open(fname.c_str(), ios_base::in);
		
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			string head_str;
			string read_seq ="";
			vector<string> tmpvec = tok.Tokenize(str,"\t" );
			map<string,string>::iterator chrmap_itr = chrmap.find(tmpvec[0]);
			
			if(chrmap_itr != chrmap.end()){
				string & chr_str = chrmap_itr -> second;
			
				int  start_pos= atoi(tmpvec[1].c_str()) ; 
				vector<string> sizevec = tok.Tokenize(tmpvec[10],"," );
				vector<string> offsetvec = tok.Tokenize(tmpvec[11],"," );
				ss.str("");
				ss<<">" << tmpvec[0] << ":";
				
				int flag =1;
				
				for(int i =0;i<sizevec.size(); i++){
					int size = atoi(sizevec[i].c_str());
					int pos = start_pos + atoi(offsetvec[i].c_str());
					
					read_seq+= chr_str.substr(pos,size);			
				 	ss<< pos+1<<",";
				 
				 	
				}
				//ss.str("");
				
				
				flag=1;
				//string headstr =">"+tmpvec[0]+":"+tmpvec[1]+":"+tmpvec[10]+":"+tmpvec[11];
				if(flag==1){
					ss<< tmpvec[3];
					string read_name = ss.str();
					
						if(tmpvec[5].compare("+") == 0){
						
							PosReadsVec.push_back(ss.str());
							readmap.insert(make_pair(tmpvec[3],""));
							//cout <<ss.str() << endl;
							PosReadsVec.push_back(read_seq);
							//cout <<read_seq << endl;
																					
						}
						else{
							//PosReadsVec.push_back(ss.str());
							
							PosReadsVec.push_back(ss.str());
							readmap.insert(make_pair(tmpvec[3],""));
							PosReadsVec.push_back(reverse_string(read_seq));
						}
						
						
						
						
						
										
				}else{
				//	cout <<"@" << endl;
				//	cout <<"@" << endl;
				}
				
			} 
		}
	}

	for(int i =0;i<PosReadsVec.size();i++){
		
		if(PosReadsVec[i][0]=='>'){
			vector<string> tokvec = tok.Tokenize(PosReadsVec[i],",");
			string read_name = tokvec[tokvec.size()-1];
			string read_name_c ;
			
			char lw = read_name[read_name.length()-1];
			
			if(lw =='1'){
				read_name_c = read_name.substr(0,read_name.length()-1)+"2";
			}
			else{
				read_name_c = read_name.substr(0,read_name.length()-1)+"1";
			}
			
			if(readmap.find(read_name_c )!= readmap.end()){
				//cout << read_name_c << " is paired with " << read_name <<endl; 
				outvec.push_back(PosReadsVec[i]);
				//cout << PosReadsVec[i]<< endl;
				outvec.push_back(PosReadsVec[i+1]);
				//cout << PosReadsVec[i+1]<< endl;
			}
		}
	}
	
	TextIO tio;
	tio.dlwrite("/rhome/ywyang/bigdata/sinica/benchmark/"+chr_dir[0]+".fa",outvec);
	filestr.close();
}

void RawReads::bed2Fasta2(string fname ,vector<string> chr_dir){
	
	map<string,string> chrmap;
	//read chr strings
	for(int i =0 ; i< chr_dir.size();i++)
	{
		//cout << "read : " << "/rhome/ywyang/bigdata/RefSeq/mm9/" +chr_dir[i]+".fa"<<endl;
		string myString = readFasta("/rhome/ywyang/bigdata/RefSeq/mm9/" +chr_dir[i]+".fa");
		//cout << "ok..." << endl;
		std::transform(myString.begin(), myString.end(), myString.begin(), ::toupper);
    chrmap.insert(make_pair(chr_dir[i], myString));
    
	}
	
//	cout << chrmap.find(chr_dir[0])->second<< endl;
	TOK tok;
	stringstream ss;
	fstream filestr;
	
	
	filestr.open(fname.c_str(), ios_base::in);
		
	string str;	 
	 while(!filestr.eof()){
	 	getline(filestr,str);
	 	if(str.length()>0){	 		
			string head_str;
			string read_seq ="";
			vector<string> tmpvec = tok.Tokenize(str,"\t" );
			map<string,string>::iterator chrmap_itr = chrmap.find(tmpvec[0]);
		
			if(chrmap_itr != chrmap.end() && tmpvec[7].compare("exon")==0){
				string & chr_str = chrmap_itr -> second;
			
				int  start_pos= atoi(tmpvec[1].c_str()) ;
				int  end_pos= atoi(tmpvec[2].c_str()) ;
				 
				
				ss.str("");
				ss<<">" << tmpvec[0] << ":"<<start_pos+1 << ","<< end_pos;
				
				read_seq+= chr_str.substr(start_pos-2,end_pos-start_pos+4);
				
				cout <<ss.str() << endl;
				cout <<read_seq << endl;
			}
		}
	}
	
	filestr.close();
	
}