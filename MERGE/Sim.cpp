#include "Sim.h"

using namespace std;

vector<string> Sim::getAllGeneList(string exprfname, string gtfdir, string idfname){
	vector<string>outVec;
  TextIO tio;
   	
//calculate RPKM


vector<vector<string> > gexprvec =  tio.dlread(exprfname,"\t");
map<string,double> exprmap;
//sum

double sum =0;

for(int i=0;i<gexprvec.size();i++){
	sum+=atof(gexprvec[i][7].c_str());
}

for(int i=0;i<gexprvec.size();i++){
	double expr = atof(gexprvec[i][7].c_str());
	double len = atof(gexprvec[i][1].c_str());
	double rpkm = (expr*1000000000)/(len*sum);
	
	exprmap.insert( make_pair(gexprvec[i][0] , rpkm ) ) ;
	//cout << gexprvec[i][0]<< ","<< rpkm<< endl; 
}
	
	
vector<vector<string> > gidvec = tio.dlread(idfname,"\t");

map<string, string> gidmap ;
map<string, double> outgmap ;

for(int i =0 ;i< gidvec.size(); i++){
	gidmap.insert(make_pair( gidvec[i][0], gidvec[i][1]));
}

for(map<string, double>::iterator expr_itr = exprmap.begin(); expr_itr != exprmap.end(); expr_itr++){
	//cout << expr_itr->first<< ","<< expr_itr->second<< endl; 	
	string gname = gidmap.find(expr_itr->first)->second;
	map<string, double>::iterator out_itr = outgmap.find(gname) ;
	if(out_itr == outgmap.end()){
		outgmap.insert(make_pair( gname, expr_itr->second ));
	}
	else{
		if(out_itr->second < expr_itr->second){
			out_itr->second = expr_itr->second;
		}
	}
	
}


map<string,double> expr_filtered_map;
for(map<string, double>::iterator g_itr = outgmap.begin(); g_itr != outgmap.end(); g_itr++){
	cout << g_itr->first<< ","<< g_itr->second<< endl; 	
	if(g_itr->second> 1.0 ){
		expr_filtered_map.insert(make_pair(g_itr->first, g_itr->second));
		
	}
}

//find genes with multiple isoforms

for(map<string, double>::iterator g_itr = expr_filtered_map.begin(); g_itr != expr_filtered_map.end(); g_itr++){
	gene g;
  cout << g_itr->first<< ","<< g_itr->second<< endl; 	
	int asm_num = g.multiGene(gtfdir+"/"+g_itr->first);	
	if(asm_num>1){
		
		outVec.push_back(g_itr->first);
	}
}

	return outVec;

}

void Sim::baseMTX(std::string lenfname ){
	
	TextIO tio;
	
	vector<vector<double> > SampleVec;
	vector<vector<string> > gtfVec = tio.dlread(lenfname, "\t");


	for(int m =0 ; m<10; m++){
		vector<double> covVec;
		//init 
		for(int i =0 ; i<gtfVec.size(); i++){
		
		int apos = atoi(gtfVec[i][0].c_str());
		int bpos = atoi(gtfVec[i][1].c_str());
		int len = bpos - apos +1;
			for(int j = 0 ; j< len ; j++){
				covVec.push_back(0.0);
			}
		}
	
		// assign 	
		int offset =0;
		int rate =10;
	
		for(int i =0 ; i<gtfVec.size(); i++){
		
		int apos = atoi(gtfVec[i][0].c_str());
		int bpos = atoi(gtfVec[i][1].c_str());
		int len = bpos - apos +1;
		/*
		if(i==gtfVec.size()-1  ){
			rate = 160;		
		}
		*/
			for(int j = 0 ; j< len ; j++){
				if( i == gtfVec.size()-1  || i==0 ){
					rate = 930;
				}else{
					rate = 930;
				}
				covVec[offset+j] = rate+ (rand() % 10)-5;
			}
			offset += len;
		}
	
		for(int i =0 ; i < covVec.size() ; i++){
		cout << covVec[i] << endl;;
		}
		//covVec[0]=0;
		SampleVec.push_back(covVec);
	}
	for(int m =0 ; m<10; m++){
		vector<double> covVec;
		//init 
		for(int i =0 ; i<gtfVec.size(); i++){
		
		int apos = atoi(gtfVec[i][0].c_str());
		int bpos = atoi(gtfVec[i][1].c_str());
		int len = bpos - apos +1;
			for(int j = 0 ; j< len ; j++){
				covVec.push_back(0.0);
			}
		}
	
		// assign 	
		int offset =0;
		int rate =10;
	
		for(int i =0 ; i<gtfVec.size(); i++){
		
		int apos = atoi(gtfVec[i][0].c_str());
		int bpos = atoi(gtfVec[i][1].c_str());
		int len = bpos - apos +1;
			for(int j = 0 ; j< len ; j++){
				if( i == gtfVec.size()-1  ){
					rate = 115;
				}else{
					rate = 115;
				}
				covVec[offset+j] = rate+ (rand() % 10)-5;
			}
				
			
			offset += len;
		}
	
		for(int i =0 ; i < covVec.size() ; i++){
		cout << covVec[i] << endl;;
		}
		//covVec[0]=0;
		SampleVec.push_back(covVec);
	}
	tio.dlwrite_double("../sim/sigfuge_data.txt",SampleVec ,"\t");
}

std::vector<std::string> Sim::simReads(std::string ExprStr, std::string BedStr ,std::string tid, int read_num , int read_len, int bias){
	//output transcript expr
	TextIO tio ; 
	TOK tok;
	stringstream ss;
	
  vector<string> outVec;
	

	vector<string> TmpExprVec ;
	TmpExprVec.push_back(ExprStr);
	tio.dlwrite( "tmpexpr.txt", TmpExprVec );

	vector<string> TmpBedVec ;
	TmpBedVec.push_back(BedStr);
	tio.dlwrite( "tmpbed.txt", TmpBedVec );
	
	
	//output reads bed
	if( bias ==1){
		
		ss.str("");
		ss << "gensimreads.py -e tmpexpr.txt -b /rhome/ywyang/bigdata/AS/sim/illuminabias.txt  -n " << read_num << " -l "<<  read_len << " -p 200,20 --stranded -o out.bed tmpbed.txt";
	}else{
		ss.str("");
		ss << "gensimreads.py -e tmpexpr.txt -n " << read_num << " -l "<<read_len <<" -p 200,20 --stranded -o out.bed tmpbed.txt";
	}
	
	string cmd = ss.str() ;
	system(cmd.c_str());
	system("touch out.bed");
	outVec =  tio.dlread_vector("out.bed");
	cmd = "rm -rf  tmpbed.txt tmpexpr.txt out.bed tmpbed.txt" ;
	system(cmd.c_str());
	//
	
	return outVec;
	
}


std::vector<std::string> Sim::simReads(std::map<std::string, std::string> ExprMap, std::map<std::string, std::string> BedMap ,std::string tid, int read_num , int len, int bias){
	//output transcript expr
	TextIO tio ; 
	TOK tok;
	stringstream ss;
	
  vector<string> outVec;
	
	if(ExprMap.find(tid)!=ExprMap.end() ){
		vector<string> TmpExprVec ;
		string exprstr = ExprMap.find(tid)->second;
		TmpExprVec.push_back(exprstr);
		tio.dlwrite( "tmpexpr.txt", TmpExprVec );
	}
	else{
		cout<< "no :" << tid <<" in ExprMap"<< endl;
		return outVec;
	}

	//output transcript bed
	if(BedMap.find(tid)!=BedMap.end() ){
		vector<string> TmpBedVec ;
		string bedstr = BedMap.find(tid)->second;
		TmpBedVec.push_back(bedstr);
		tio.dlwrite( "tmpbed.txt", TmpBedVec );
	}else{
		cout<< "no :" << tid <<" in BedMap"<< endl;
		return outVec;
	}
	
	//output reads bed
	if( bias ==1){
		
		ss.str("");
		ss << "gensimreads.py -e tmpexpr.txt -b /rhome/ywyang/bigdata/AS/sim/illuminabias.txt  -n " << read_num << " -l "<<len << " -p 200,20 --stranded -o out.bed tmpbed.txt";
	}else{
		ss.str("");
		ss << "gensimreads.py -e tmpexpr.txt -n " << read_num << " -l "<< len<<" -p 200,20 --stranded -o out.bed tmpbed.txt";
	}
	
	string cmd = ss.str() ;
	system(cmd.c_str());
	system("touch out.bed");
	outVec =  tio.dlread_vector("out.bed");
	cmd = "rm -rf  tmpbed.txt tmpexpr.txt out.bed tmpbed.txt" ;
	system(cmd.c_str());
	//
	
	return outVec;
	
}

std::vector<std::string> Sim::simReads_para(std::map<std::string, std::string> ExprMap, std::map<std::string, std::string> BedMap ,std::string tid, int read_num , int len, int bias, int IDX){
	//output transcript expr
	TextIO tio ; 
	TOK tok;
	stringstream ss;
	
  vector<string> outVec;
	
	if(ExprMap.find(tid)!=ExprMap.end() ){
		vector<string> TmpExprVec ;
		string exprstr = ExprMap.find(tid)->second;
		TmpExprVec.push_back(exprstr);
		ss.str("");
		ss<< "tmpexpr.txt"<< IDX;
		tio.dlwrite( ss.str(), TmpExprVec );
	}
	else{
		cout<< "no :" << tid <<" in ExprMap"<< endl;
		return outVec;
	}

	//output transcript bed
	if(BedMap.find(tid)!=BedMap.end() ){
		vector<string> TmpBedVec ;
		string bedstr = BedMap.find(tid)->second;
		TmpBedVec.push_back(bedstr);
		ss.str("");
		ss<< "tmpbed.txt" << IDX;
		tio.dlwrite( ss.str(), TmpBedVec );
	}else{
		cout<< "no :" << tid <<" in BedMap"<< endl;
		return outVec;
	}
	
	//output reads bed
	if( bias ==1){
		
		ss.str("");
		ss << "python /rhome/ywyang/bigdata/3rd/RNASeqReadSimulator/src/gensimreads.py -e tmpexpr.txt"<< IDX<<" -b /rhome/ywyang/bigdata/AS/sim/illuminabias.txt  -n " << read_num << " -l "<<len << " -p 200,20 --stranded -o out.bed"<< IDX<<" tmpbed.txt"<< IDX;
	}else{
		ss.str("");
		ss << "python /rhome/ywyang/bigdata/3rd/RNASeqReadSimulator/src/gensimreads.py -e tmpexpr.txt"<< IDX<< " -n " << read_num << " -l "<< len<<" -p 200,20 --stranded -o out.bed"<< IDX<< " tmpbed.txt"<<IDX;
	}
	
	string cmd = ss.str() ;
	system(cmd.c_str());
	
	ss.str("");
	ss<<"touch out.bed"<<IDX ;
	cmd = ss.str();
	system(cmd.c_str());
	
	ss.str("");
	ss<<"out.bed"<<IDX;
	outVec =  tio.dlread_vector(ss.str());
	
	ss.str("");
	ss<< "rm -rf "<<"tmpbed.txt"<< IDX<< " tmpexpr.txt"<< IDX<< " out.bed"<< IDX<< " tmpbed.txt"<< IDX ;
	cmd=ss.str();
	system(cmd.c_str());
	//
	
	return outVec;
	
}

std::map< std::string , double > Sim::calc_RC(std::string exprfname){
	map< string , double > outMap ; 
	TextIO tio ; 
	vector<vector<string> > datavec = tio.dlread(exprfname,  "\t" );
	
	
	//calc all weighted length
	double all_value = 0;
	
	for(int i =0 ; i< datavec.size() ;i++){
		
		double weight = atof(datavec[i][7].c_str());
		all_value += weight;
	}
	
	//cout << all_value << endl;
	all_value =   all_value/1000000;
	for(int i =0 ; i< datavec.size() ;i++){
		
		double weight = atof(datavec[i][7].c_str());
		outMap.insert(make_pair(datavec[i][0],weight/all_value));
	}

	
	return outMap ; 	
}


std::map< std::string , double > Sim::calc_RC_sp(){
	
	map< string , double > outMap ; 
	map< string , double > ExprMap = calc_RC("/rhome/ywyang/bigdata/AS/sim/expr/1.expr");
	map<string , string> g_tid_map =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	TOK tok ;
	TextIO tio;
	vector<string> spvec = tio.dlread_vector( "/rhome/ywyang/bigdata/AS/sim/gene_list/sp_genes" ); 
	
	map<string, string> spmap ; 
	for(int i =0 ; i< spvec.size() ; i++){
		spmap.insert(make_pair(spvec[i] , ""));
	}
	
	vector<string> allvec = tio.dlread_vector( "/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt" ); 
	
	map<string, string> allmap ; 
	for(int i =0 ; i< allvec.size() ; i++){
		allmap.insert(make_pair(allvec[i] , ""));
	}

	
	map<string,string> usedmap;
	
	vector<string> idvec  = tio.dlread_vector("/rhome/ywyang/bigdata/RefSeq/hg38/Transcript_Gene_ID");
	map<string , string> tid_g_map ; 
	for(int i =0; i<idvec.size(); i++ ){
		vector<string> tokvec = tok.Tokenize(idvec[i], "\t");
		tid_g_map.insert(make_pair( tokvec[0],tokvec[1] ) ); 
	}
		
	for( map< string , double >::iterator d_itr = ExprMap.begin() ; d_itr !=ExprMap.end() ;d_itr++ ){
	
		string tid = d_itr->first;
		string gname = tid_g_map.find(tid)->second;
		if(allmap.find(gname) != allmap.end() ){
			if(spmap.find(gname)  != spmap.end()  ){
				if(usedmap.find(gname)==usedmap.end()){
			//print all transcripts of the gene
			
				vector<string> tokvec = tok.Tokenize(g_tid_map.find(gname)->second,"_");
			
				GenomeSegment gs ;
				gs.initGenomeSegment("/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/"+gname,INT_MAX);
			
				for(int j=1; j< tokvec.size() ; j++){
				
				
				int g_len = 0;
							
				for( map<string , map<string,string> >::iterator g_itr = gs.ExonMap.begin() ;g_itr != gs.ExonMap.end() ; g_itr ++){
					if((g_itr->second).find("\""+tokvec[j]+"\"") != (g_itr->second).end()){
						vector<string> tokvec = tok.Tokenize(g_itr->first, "_");
						int len = atoi(tokvec[1].c_str())-atoi(tokvec[0].c_str())+1; 
									//cout << g_itr->first<< ","<< len<< endl;
						g_len+= len;
					}
				}
				//cout <<g_len<< ":";
				
				//cout << tokvec[j]<< "\t";
				double mean =0;
		
				if(j==1){
							mean = 40.0*((double)g_len/(double)1000);
								
				}
				if(j==2){
						mean = 3.0*((double)g_len/(double)1000);
				}
				if(j>2){
						mean = 0.3*((double)g_len/(double)1000);
				}
				
					
					if(mean > 10000){
							mean = 10000;
					}
					cout << tokvec[j]<< "\t" << mean<< endl;
				}
			
			//cout << endl;
			
			usedmap.insert(make_pair(gname,""));
				}
			}
			else{
				
			//cout <<gname << endl;
			if(40 *d_itr->second < 10000){
				cout << d_itr->first << "\t"<< 40 *d_itr->second <<endl;
			}
			else{
				cout << d_itr->first << "\t"<< 10000 <<endl;
			}
			outMap.insert( make_pair(d_itr->first ,d_itr->second) );
			}
		}
	}	
	
	return outMap ; 	
}

std::map< std::string , std::string > Sim::read_transcripts_map(std::string glist_fname , std::string gtfpath ){
	TextIO tio;
	map< string,string > outMap;
	GenomeSegment gs;
	stringstream ss;
		
	vector<vector<string> > gidvec = tio.dlread(glist_fname,"\t");
	
	for(int i =0 ; i < gidvec.size() ; i++){
		string gname = gidvec[i][0];
		
		vector<vector<string> > gtfvec = tio.dlread(gtfpath+"/"+gname,"\t");
		string chrstr;
		map<string, string> tidmap;
		
		for(int j =0 ; j<gtfvec.size() ; j++){
			chrstr = gtfvec[j][0];
			map<string,string> tmpmap = gs.parse_gid(gtfvec[j][8]);
			string tid = tmpmap.find("transcript_id")->second;
			if(tidmap.find(tid)== tidmap.end()){
				
				tidmap.insert(make_pair( tid.substr(1,tid.length()-2) , "" ));
			}
		}
		
		map<string, string>::iterator mapitr;
		ss.str("");
		ss<< chrstr;
		for(mapitr  = tidmap.begin() ;mapitr  != tidmap.end() ; mapitr++){
				 ss<< "_"<< mapitr->first ;
		}
		outMap.insert(make_pair( gname, ss.str() ));
					
	}
	
	
	
	return outMap;
}

void Sim::update_RC(){
	TextIO tio;
	TOK tok;
	map<string , string> tidmap =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	
	vector<string> gvec = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/gene_list/DE_genes_new");
	vector<string> allgvec = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt");
	
	map<string,double> rcmap1,rcmap2; 
	map<string,double> dmap,scalemap;
	//init rcmap1 
	
	vector<string> rctmpvec = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/expr/binary/1.rc");
	for(int i =0; i<rctmpvec.size(); i++){
		vector<string> tmpvec = tok.Tokenize( rctmpvec[i], "\t");
		rcmap1.insert(make_pair(tmpvec[0],atof(tmpvec[1].c_str()) ));
	}
	
	//init rcmap2
	
	rctmpvec = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/expr/binary/2.rc");
	for(int i =0; i<rctmpvec.size(); i++){
		vector<string> tmpvec = tok.Tokenize( rctmpvec[i], "\t");
		rcmap2.insert(make_pair(tmpvec[0],atof(tmpvec[1].c_str()) ));
	}
	
	//init gene_pred
	map<string,string> predMap;
	rctmpvec = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/sim_reads/pred/sdeap_20_4/genes_pred.txt");
	for(int i =0; i<rctmpvec.size(); i++){
		vector<string> tmpvec = tok.Tokenize( rctmpvec[i], "\t");
		predMap.insert(make_pair(tmpvec[0], rctmpvec[i] ));
	}
	 
	for(int i=0 ; i<gvec.size(); i++){
		string gname = gvec[i];
		cout << gname<<"\t" ;
		string exstr = predMap[gname];
		vector<string> tmpvec = tok.Tokenize( exstr, "\t");
		int size = (tmpvec.size() -5)/2;
		double sum =0;
		for(int j = 0; j<size ; j++){
			double d = atof(tmpvec[j+3].c_str()); 
			//cout << d<<",";
			sum+= d; 
		}
		//cout << endl;
		double avg = sum/size;
		double scale = (3*((double) rand() / (RAND_MAX))+5.5)/avg ;
		if(avg > 0){
			
			
			if(scale > 1.0 ){
				
			}
			else{
				scale =1;
			}
			
		}
		else{
			scale = 1;
		}
		scalemap.insert(make_pair( gname, scale));
		//cout<< avg*scale<< ","<< avg << ","<< scale << endl;
		cout << avg << "\t" << scale << endl;
	}
	/*
	stringstream ss;
	vector<string> r1vec,r2vec ; 
	
	
	for(int i =0; i<allgvec.size(); i++){
		string gname = allgvec[i];
		map<string,string>::iterator name_itr = tidmap.find(gname);
		
		vector<string> tokvec = tok.Tokenize(name_itr->second , "_");
		double scale =1;
		if(scalemap.find(gname)!=scalemap.end()){
			scale =scalemap[gname];
			cout<< gname << endl;
			cout<< name_itr->second << endl;
		}
		string chr_name = tokvec[0];
		
		for(int j=1; j<tokvec.size(); j++){
			ss.str("");
			ss<< tokvec[j]<< "\t"<< scale*rcmap1[tokvec[j]];
			
			r1vec.push_back(ss.str()); 		 
			
			ss.str("");
			ss<< tokvec[j]<< "\t"<< scale*rcmap2[tokvec[j]];
			r2vec.push_back(ss.str()); 		 
		
		}
	}	
	tio.dlwrite( "/rhome/ywyang/bigdata/AS/sim/expr/binary/1_new.rc" , r1vec);
	tio.dlwrite( "/rhome/ywyang/bigdata/AS/sim/expr/binary/2_new.rc" , r2vec);
	*/
}

void Sim::simSamples_para(std::vector<int>  GroupVec , std::string rcfname	, double lib_size ,double dispersion ,	std::string reads_outpath , int IDX ){
	
	stringstream ss;
	
	map<string, double>::iterator rc_itr ; 	
	
	//read gene-transcript map
	cout << "read gene-transcript map"<< endl; 
  map<string , string> tidmap =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	map<string , string>::iterator name_itr ;
	TextIO tio;
	TOK tok;
	//vector<vector<string> > testgenesvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt","\t");
	vector<vector<string> > testgenesvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/test_genes","\t");
	//read transcript-bed map by chrs;
	cout << "read transcript-bed map by chrs" << endl;
	map<string, map<string,string> > ChrBedMap;
	for(int i =1;i<=23 ; i++){
		string chrname ;
		if(i==23){
			chrname = "chrX";
			
		}else{
			ss.str("");
			ss<<"chr"<<i;
			chrname =ss.str();
		}
		vector<string> bedvec = tio.dlread_vector("/rhome/ywyang/bigdata/RefSeq/hg38/chr_bed/"+chrname);
		map<string,string> bedmap;
		for(int j =0 ; j<bedvec.size(); j++){
			 
			vector<string> tokvec = tok.Tokenize( bedvec[j], "\t");
			string tid = tokvec[3];
			bedmap.insert(make_pair( tid , bedvec[j] ));
			
		}
		
		ChrBedMap.insert(make_pair(chrname,bedmap));
	}
	//read ExprMap vector
	cout << "read ExprMap vector" << endl;
	vector<map<string,string> > ExprMapVec;
	for(int n=0;n<GroupVec.size();n++ ){
		map<string,string> ExprMap;
		ss.str("");
		ss<< "/rhome/ywyang/bigdata/AS/sim/expr/1.expr" ;
		
		vector<string> bedvec = tio.dlread_vector(ss.str());
		for(int i =0 ; i<bedvec.size(); i++){
			vector<string> tokvec  = tok.Tokenize(bedvec[i] ,"\t" );
			string tid = tokvec[0];
			ExprMap.insert(make_pair(tid , bedvec[i] ));
		}
		ExprMapVec.push_back(ExprMap);
	}
	
	//read RCMap 
	cout << "read RCMap vector" << endl;
	map<string , vector<string> > rcmap;
	vector<string> rctmpvec = tio.dlread_vector(rcfname);
	for(int i =0; i<rctmpvec.size(); i++){
		vector<string> tmpvec = tok.Tokenize( rctmpvec[i], "\t");
		rcmap.insert(make_pair(tmpvec[0],tmpvec));
	}
	
	
	//
	int offset =0; 
	//generate reads for m samples of n groups
	for(int n=0 ;n<GroupVec.size() ; n++){
		int sample_num = GroupVec[n];
		ss.str("");
		if(n>0){
			offset += GroupVec[n-1];
		}
	//	ss<< "/rhome/ywyang/bigdata/AS/sim/expr/"<< n+1<<".expr" ;
	//	map<string, double> rcmap =  calc_RC(ss.str());
		
		for(int m =0 ; m<sample_num; m++){
			for(int i =0 ;  i < testgenesvec.size() ; i++){
				//if(i%30 == IDX){
				
				vector<string> totalreads;
				string gname = testgenesvec[i][0];
				cout << "process...."<< gname <<endl;
				name_itr = tidmap.find(gname);
				if(name_itr != tidmap.end()){
					vector<string> tokvec = tok.Tokenize(name_itr->second , "_");
					string chr_name = tokvec[0];
				//for every transcript
					for(int j=1;j<tokvec.size() ; j++ ){
						string tid = tokvec[j];
						vector<string> rcvec = rcmap.find(tid)->second;
					//	rc_itr = rcmap.find(tid);
					//	cout << rc_itr->first <<","<<  rc_itr-> second*lib_size << endl;
						double mean  = atoi(rcvec[offset+m+1].c_str());
						//cout <<"mean: " << mean << endl;
						
						
						//access rcnum ;
						int rc_num = mean;
													
						cout << tid<< ", " << mean<< endl;
						map<string, string> BedMap = (ChrBedMap.find(chr_name))->second;
							
						vector<string> readvec =simReads_para(ExprMapVec[n], BedMap,tid,   rc_num, 50, 1,IDX);
						for(int k =0 ; k<readvec.size(); k++){
							totalreads.push_back(readvec[k]);
						}	
							
					
					
					}
					ss.str("");
					ss<< reads_outpath<< "/"<< gname<< "_" << n<< "_" <<m<<".bed";  
					tio.dlwrite(ss.str() , totalreads);
				
				}
		
				//}
			}	
		}
		
	}

}

void Sim::simSamples(std::vector<int>  GroupVec ,	double lib_size ,double dispersion ,	std::string reads_outpath ){
	
	stringstream ss;
	
	map<string, double>::iterator rc_itr ; 	
	
	//read gene-transcript map
//	map<string , string> tidmap =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	map<string , string> tidmap =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/test_genes", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	map<string , string>::iterator name_itr ;
	TextIO tio;
	TOK tok;
//	vector<vector<string> > testgenesvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt","\t");
	vector<vector<string> > testgenesvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/test_genes","\t");
	//read transcript-bed map by chrs;
	
	map<string, map<string,string> > ChrBedMap;
	for(int i =1;i<=23 ; i++){
		string chrname ;
		if(i==23){
			chrname = "chrX";
			
		}else{
			ss.str("");
			ss<<"chr"<<i;
			chrname =ss.str();
		}
		vector<string> bedvec = tio.dlread_vector("/rhome/ywyang/bigdata/RefSeq/hg38/chr_bed/"+chrname);
		map<string,string> bedmap;
		for(int j =0 ; j<bedvec.size(); j++){
			 
			vector<string> tokvec = tok.Tokenize( bedvec[j], "\t");
			string tid = tokvec[3];
			bedmap.insert(make_pair( tid , bedvec[j] ));
			
		}
		
		ChrBedMap.insert(make_pair(chrname,bedmap));
	}
	//read ExprMap vector
	vector<map<string,string> > ExprMapVec;
	for(int n=0;n<GroupVec.size();n++ ){
		map<string,string> ExprMap;
		ss.str("");
		ss<< "/rhome/ywyang/bigdata/AS/sim/expr/"<< n+1<<".expr" ;
		
		vector<string> bedvec = tio.dlread_vector(ss.str());
		for(int i =0 ; i<bedvec.size(); i++){
			vector<string> tokvec  = tok.Tokenize(bedvec[i] ,"\t" );
			string tid = tokvec[0];
			ExprMap.insert(make_pair(tid , bedvec[i] ));
		}
		ExprMapVec.push_back(ExprMap);
	}
	
	//read RCMap 
	
	//generate reads for m samples of n groups
	for(int n=0 ;n<GroupVec.size() ; n++){
		int sample_num = GroupVec[n];
		ss.str("");
		ss<< "/rhome/ywyang/bigdata/AS/sim/expr/"<< n+1<<".expr" ;
		map<string, double> rcmap =  calc_RC(ss.str());
		
		for(int m =0 ; m<sample_num; m++){
			for(int i =0 ;  i < testgenesvec.size() ; i++){
			vector<string> totalreads;
			string gname = testgenesvec[i][0];
			
			name_itr = tidmap.find(gname);
				if(name_itr != tidmap.end()){
					vector<string> tokvec = tok.Tokenize(name_itr->second , "_");
					string chr_name = tokvec[0];
				//for every transcript
					for(int j=1;j<tokvec.size() ; j++ ){
						string tid = tokvec[j];
						rc_itr = rcmap.find(tid);
						cout << rc_itr->first <<","<<  rc_itr-> second*lib_size << endl;
						double mean  = rc_itr-> second*lib_size;
						cout <<"mean: " << mean << endl;
						if(mean > 10000){
							mean = 10000;
						}
						if(mean > 1.0){
						
							//usr RCGen.r tp generate the read count
							ss.str("");
							ss<< "Rscript RCGen.r " << mean << " " << dispersion;
							string cmd = ss.str();
							system(cmd.c_str());
							vector<string> rcvec = tio.dlread_vector("rctmp");
							int rc_num = atoi(rcvec[0].c_str());							
							
							map<string, string> BedMap = (ChrBedMap.find(chr_name))->second;
							
							vector<string> readvec =simReads(ExprMapVec[n], BedMap,tid,  rc_num, 50, 1);
							for(int k =0 ; k<readvec.size(); k++){
								totalreads.push_back(readvec[k]);
							}	
							
						}
					
					}
					ss.str("");
					ss<< reads_outpath<< "/"<< gname<< "_" << n<< "_" <<m<<".bed";  
					tio.dlwrite(ss.str() , totalreads);
				
				}
		
			}	
		}
		
	}

}

void Sim::simSamples_test(std::vector<int>  GroupVec ,	double lib_size ,double dispersion ,	std::string reads_outpath ){
	
	stringstream ss;
	
	map<string, double>::iterator rc_itr ; 	
	
	//read gene-transcript map
	//map<string , string> tidmap =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	map<string , string> tidmap =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/sp_genes", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	map<string , string>::iterator name_itr ;
	TextIO tio;
	TOK tok;
	vector<vector<string> > testgenesvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/sp_genes","\t");
	
	//read transcript-bed map by chrs;
	
	map<string, map<string,string> > ChrBedMap;
	for(int i =1;i<=23 ; i++){
		string chrname ;
		if(i==23){
			chrname = "chrX";
			
		}else{
			ss.str("");
			ss<<"chr"<<i;
			chrname =ss.str();
		}
		vector<string> bedvec = tio.dlread_vector("/rhome/ywyang/bigdata/RefSeq/hg38/chr_bed/"+chrname);
		map<string,string> bedmap;
		for(int j =0 ; j<bedvec.size(); j++){
			 
			vector<string> tokvec = tok.Tokenize( bedvec[j], "\t");
			string tid = tokvec[3];
			bedmap.insert(make_pair( tid , bedvec[j] ));
			
		}
		
		ChrBedMap.insert(make_pair(chrname,bedmap));
	}
	
	vector<map<string,string> > ExprMapVec;
		for(int n=0;n<GroupVec.size();n++ ){
		map<string,string> ExprMap;
		ss.str("");
		ss<< "/rhome/ywyang/bigdata/AS/sim/expr/"<< n+1<<".expr" ;
		
		vector<string> bedvec = tio.dlread_vector(ss.str());
		for(int i =0 ; i<bedvec.size(); i++){
			vector<string> tokvec  = tok.Tokenize(bedvec[i] ,"\t" );
			string tid = tokvec[0];
			ExprMap.insert(make_pair(tid , bedvec[i] ));
		}
		ExprMapVec.push_back(ExprMap);
	}
	//read RCMap 
	
	//generate reads for m samples of n groups
	for(int n=0 ;n<GroupVec.size() ; n++){
		int sample_num = GroupVec[n];
		ss.str("");
		ss<< "/rhome/ywyang/bigdata/AS/sim/expr/"<< n+1<<".expr" ;
		map<string, double> rcmap =  calc_RC(ss.str());
		
		for(int m =0 ; m<sample_num; m++){
			for(int i =0 ;  i <testgenesvec.size() ; i++){
			vector<string> totalreads;
			string gname = testgenesvec[i][0];
			
			name_itr = tidmap.find(gname);
				if(name_itr != tidmap.end()){
					vector<string> tokvec = tok.Tokenize(name_itr->second , "_");
					string chr_name = tokvec[0];
				//for every transcript
					for(int j=1;j<2/*tokvec.size()*/ ; j++ ){
						
						string tid = tokvec[j];
						rc_itr = rcmap.find(tid);
						cout << rc_itr->first <<","<<  rc_itr-> second*lib_size << endl;
						double mean =0;
						if(i<0){
						  mean = rc_itr-> second*lib_size;
						}
						else{
							//calc gene length 
							GenomeSegment gs ;
							gs.initGenomeSegment("/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/"+gname,INT_MAX);
							//gs._debug();
							int g_len = 0;
							
							for( map<string , map<string,string> >::iterator g_itr = gs.ExonMap.begin() ;g_itr != gs.ExonMap.end() ; g_itr ++){
								if((g_itr->second).size() >0){
									vector<string> tokvec = tok.Tokenize(g_itr->first, "_");
									
									int len = atoi(tokvec[1].c_str())-atoi(tokvec[0].c_str())+1; 
									//cout << g_itr->first<< ","<< len<< endl;
									g_len+= len;
								}
							}
							
							//
							
							if(n==0){
								mean = 40.0*((double)g_len/(double)1000);
								
							}
							else{
								mean = 11.0*((double)g_len/(double)1000);
							}
						}
						cout <<"mean: " << mean << endl;
						int rc_num =0;
						if(mean > 1.0 ){
						
							//usr RCGen.r tp generate the read count
							ss.str("");
							ss<< "Rscript RCGen.r " << mean << " " << dispersion;
							string cmd = ss.str();
							
							system(cmd.c_str());
							vector<string> rcvec = tio.dlread_vector("rctmp");
							rc_num = atoi(rcvec[0].c_str());							
							
							
						}
						else{
							rc_num =1;
						}
						
						cout << chr_name << endl;
						map<string, string> BedMap = (ChrBedMap.find(chr_name))->second;	
						
						vector<string> readvec =simReads(ExprMapVec[n], BedMap,tid,  rc_num, 50, 1);
						for(int k =0 ; k<readvec.size(); k++){
							totalreads.push_back(readvec[k]);
						}	
							
							
						
					
					}
					ss.str("");
					ss<< reads_outpath<< "/"<< gname<< "_" << n<< "_" <<m<<".bed";  
					tio.dlwrite(ss.str() , totalreads);
				
				}
		
			}	
		}
		
	}

}



std::vector<std::string> Sim::createDSgenes(std::string exprfname, std::vector<std::string> testgeneVec, std::string idfname){
  vector<string> outVec;
	TextIO tio;
	TOK tok;
	stringstream ss;
	//load gid
	map<string, string> gidmap  ; 
	vector<vector<string> > gidvec = tio.dlread(idfname,"\t");
	
	for(int i =0 ;i< gidvec.size(); i++){
		gidmap.insert(make_pair( gidvec[i][0], gidvec[i][1]));
	}
	
	//load the expr file

	vector<vector<string> > gexprvec =  tio.dlread(exprfname,"\t");
	
	map<string,double> exprmap;
	map<string,vector<string> > gene_RNA_exprmap;
	//sum

	double sum =0;

	for(int i=0;i<gexprvec.size();i++){
		sum+=atof(gexprvec[i][7].c_str());
	}

	for(int i=0;i<gexprvec.size();i++){
		
		double expr = atof(gexprvec[i][7].c_str());
		double len = atof(gexprvec[i][1].c_str());
		double rpkm = (expr*1000000000)/(len*sum);
	
		exprmap.insert( make_pair(gexprvec[i][0] , rpkm ) ) ;
		map<string , string>::iterator gid_itr = gidmap.find(gexprvec[i][0]);
		string gname = 	gid_itr->second;
		map<string,vector<string> >::iterator gene_RNA_itr = gene_RNA_exprmap.find(gname);
		
		if(gene_RNA_itr == gene_RNA_exprmap.end()){
			vector<string> tmpvec ; 
			ss.str("");
			ss<< gexprvec[i][0]<< ","<< rpkm;
			tmpvec.push_back(ss.str());
			gene_RNA_exprmap.insert(make_pair(gname, tmpvec));
					
		}
		else{
			ss.str("");
			ss<< gexprvec[i][0]<< ","<< rpkm;
			(gene_RNA_itr->second).push_back(ss.str());
			
		}
	//cout << gexprvec[i][0]<< ","<< rpkm<< endl; 
	}
	

  map<string,string > ShuffleMap ; 
	//find gene with isoforms expressed significantly (> 1/4)
	
	
	for(int i = 0; i< testgeneVec.size() ; i++){
		
		map<string,vector<string> >::iterator expr_itr = gene_RNA_exprmap.find(testgeneVec[i]);
		vector<string> 	isovec = expr_itr->second;
		double min =  std::numeric_limits<double>::max();
		int min_idx = -1;
		double max = std::numeric_limits<double>::min();
		int max_idx = -1;	
		
		for(int j =0; j<isovec.size();j++){
			vector<string> tokvec = tok.Tokenize(isovec[j],",");
			double rpkm = atof(tokvec[1].c_str()); 
			if(rpkm>max){
				max = rpkm;
				max_idx = j;
			}
			
			if(rpkm<min){
				min = rpkm;
				min_idx = j;
			}
			//cout << testgeneVec[i]<< ":"<< isovec[j]<<endl;
		}
		vector<string> tokvec = tok.Tokenize(isovec[max_idx],",");
		double max_rpkm =  atof(tokvec[1].c_str());
		tokvec = tok.Tokenize(isovec[min_idx],",");
		double min_rpkm = atof(tokvec[1].c_str());
		
		if((max_rpkm/min_rpkm)>10.0){
			
			ShuffleMap.insert(make_pair( testgeneVec[i],isovec[max_idx] +","+ isovec[min_idx]));
		}
		else{
			//cout <<testgeneVec[i]<< ","<< isovec[max_idx]<< ","<< max_rpkm << ","<< isovec[min_idx]<< "," << min_rpkm<< endl;
		}
	}

	
	
	
	int i=0;
	//select 10% genes
	
	vector<string> dsvec;
	
	for( map<string,string >::iterator s_itr = ShuffleMap.begin() ;s_itr != ShuffleMap.end() ; s_itr++ ){
		int flag = 0;
		if(i%10==0 || i<4){
			dsvec.push_back(s_itr->first);
			flag = 1;
		}
		i++;
		
	}

	
	//
  map<string,string> upmap,downmap,dsmap;
  vector<string> upvec , downvec, svec;
  for(i =0 ; i<dsvec.size() ; i++){
  	
  	if(i%3 == 0){
  		vector<string> tokvec = tok.Tokenize(ShuffleMap.find(dsvec[i])->second,",");
  		upvec.push_back(dsvec[i]);
  		upmap.insert(make_pair(tokvec[0],""));
  		
  		
  	}
  	if(i%3 == 1){
  		vector<string> tokvec = tok.Tokenize(ShuffleMap.find(dsvec[i])->second,",");
  		downvec.push_back(dsvec[i]);
  		downmap.insert(make_pair(tokvec[0],""));
  		
  	}
  	if(i%3 == 2){
  		vector<string> tokvec = tok.Tokenize(ShuffleMap.find(dsvec[i])->second,",");
  		svec.push_back(dsvec[i]);
  		dsmap.insert(make_pair(tokvec[0],tokvec[2]));
  		dsmap.insert(make_pair(tokvec[2],tokvec[0]));
  		
  	}
  	
  }
  
  tio.dlwrite("up_genes", upvec);
	tio.dlwrite("down_genes", downvec);
	tio.dlwrite("ds_genes", svec);
	//expr 
	/*
	map<string, string> exmap ;  
	for(int i =0 ; i <gexprvec.size() ; i++){
		exmap.insert(make_pair( gexprvec[i][0], gexprvec[i][7]));
	}
	ss.precision(15);
	for(int i =0 ; i <gexprvec.size() ; i++){
		
		if(dsmap.find(gexprvec[i][0])!= dsmap.end()){
			string bname = dsmap.find(gexprvec[i][0])->second;
			//cout << gexprvec[i][0]<< "<->"<< bname<< endl;
			string b_rpkm = exmap.find(bname)->second; 
			ss.str("");
			ss << gexprvec[i][0]  << "\t"<< gexprvec[i][1]<< "\t"<< gexprvec[i][2]<< "\t"<< gexprvec[i][3]<< "\t"<<gexprvec[i][4] << "\t" << gexprvec[i][5] << "\t"<< gexprvec[i][6]<< "\t"<< b_rpkm ;
			outVec.push_back(ss.str());
		} else if(upmap.find(gexprvec[i][0])!= upmap.end()){
			double r = 3*((double) rand() / (RAND_MAX))+2;
			double rate = pow(2,r); 
			ss.str("");
			ss << gexprvec[i][0] << "\t"<< gexprvec[i][1]<< "\t"<< gexprvec[i][2]<< "\t"<< gexprvec[i][3]<< "\t"<<gexprvec[i][4] << "\t" << gexprvec[i][5] << "\t"<< gexprvec[i][6]<< "\t"<< atof(gexprvec[i][7].c_str())*rate ;
			outVec.push_back(ss.str());
		}
		else if (downmap.find(gexprvec[i][0])!= downmap.end()){
			ss.str("");
			double r = 1*((double) rand() / (RAND_MAX))+4;
			double rate = pow(2,r);
			ss << gexprvec[i][0] << "\t"<< gexprvec[i][1]<< "\t"<< gexprvec[i][2]<< "\t"<< gexprvec[i][3]<< "\t"<<gexprvec[i][4] << "\t" << gexprvec[i][5] << "\t"<< gexprvec[i][6]<< "\t"<< atof(gexprvec[i][7].c_str())/rate ;
			outVec.push_back(ss.str());
		}
		else{
			ss.str("");
			ss << gexprvec[i][0] << "\t"<< gexprvec[i][1]<< "\t"<< gexprvec[i][2]<< "\t"<< gexprvec[i][3]<< "\t"<<gexprvec[i][4] << "\t" << gexprvec[i][5] << "\t"<< gexprvec[i][6]<< "\t"<< gexprvec[i][7] ;
			outVec.push_back(ss.str());
		}
		
	}
	*/
	return outVec;
}

std::vector<std::string> Sim::binary_EE_RC(std::string glistfname){
	vector<string> outVec;
	
	TextIO tio;
	TOK tok;
	stringstream ss;
	ss.precision(15);
	map<string , string> tid_g_map ;
	map<string , string> g_tid_map =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	map<string, double> exprmap;
	vector<string> idvec  = tio.dlread_vector("/rhome/ywyang/bigdata/RefSeq/hg38/Transcript_Gene_ID");
	for(int i =0; i<idvec.size(); i++ ){
		vector<string> tokvec = tok.Tokenize(idvec[i], "\t");
		tid_g_map.insert(make_pair( tokvec[0],tokvec[1] ) ); 
	}
	
	vector<string> exprvec  = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/expr/binary/1.rc");
	for(int i =0; i<exprvec.size(); i++ ){
		vector<string> tokvec = tok.Tokenize(exprvec[i], "\t");
		exprmap.insert(make_pair( tokvec[0],atof(tokvec[1].c_str())) ); 
	}
	
	
	vector<string> glistvec  = tio.dlread_vector(glistfname);
	for(int i =0; i<glistvec.size(); i++ ){
		
		vector<string> tokvec = tok.Tokenize( g_tid_map.find(glistvec[i])->second, "_");
		
		//for every transcript 
		for(int j =1 ;j<tokvec.size() ; j++){
			string tid = tokvec[j];
			ss.str("");
			ss << tid<<"\t" << exprmap.find(tid)->second ;
			outVec.push_back(ss.str());
		}
		
	}
	
	
	return  outVec;
}

std::vector<std::string> Sim::binary_up_RC(std::string glistfname,std::string spfname ){
	vector<string> outVec;
	
	TextIO tio;
	TOK tok;
	stringstream ss;
	ss.precision(15);
	map<string , string> tid_g_map ;
	map<string , string> g_tid_map =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	map<string, double> exprmap;
	vector<string> idvec  = tio.dlread_vector("/rhome/ywyang/bigdata/RefSeq/hg38/Transcript_Gene_ID");
	for(int i =0; i<idvec.size(); i++ ){
		vector<string> tokvec = tok.Tokenize(idvec[i], "\t");
		tid_g_map.insert(make_pair( tokvec[0],tokvec[1] ) ); 
	}
	
	vector<string> exprvec  = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/expr/binary/1.rc");
	for(int i =0; i<exprvec.size(); i++ ){
		vector<string> tokvec = tok.Tokenize(exprvec[i], "\t");
		exprmap.insert(make_pair( tokvec[0],atof(tokvec[1].c_str())) ); 
	}
	
	map<string, double> spmap;
	vector<string> spvec  = tio.dlread_vector(spfname);
	for(int i =0; i<spvec.size(); i++ ){
		
		spmap.insert(make_pair( spvec[0],0) ); 
	}
	
	
	vector<string> glistvec  = tio.dlread_vector(glistfname);
	for(int i =0; i<glistvec.size(); i++ ){
		
		vector<string> tokvec = tok.Tokenize( g_tid_map.find(glistvec[i])->second, "_");
		
		//for every transcript 
		for(int j =1 ;j<tokvec.size() ; j++){
			string tid = tokvec[j];
			if(spmap.find(glistvec[i]) == spmap.end()){
				ss.str("");
				double r = 3*((double) rand() / (RAND_MAX))+2;
				double rate = pow(2,r);
				if(rate > 5){
					rate =5;
				}
				ss << tid<<"\t" << rate*exprmap.find(tid)->second ;
			}
			else{
				ss.str("");
				ss << tid<<"\t" << 4*exprmap.find(tid)->second ;
			}
			outVec.push_back(ss.str());
		}
		
	}
	
	
	return  outVec;
}

std::vector<std::string> Sim::binary_down_RC(std::string glistfname,std::string spfname ){
	vector<string> outVec;
	
	TextIO tio;
	TOK tok;
	stringstream ss;
	ss.precision(15);
	map<string , string> tid_g_map ;
	map<string , string> g_tid_map =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	map<string, double> exprmap;
	vector<string> idvec  = tio.dlread_vector("/rhome/ywyang/bigdata/RefSeq/hg38/Transcript_Gene_ID");
	for(int i =0; i<idvec.size(); i++ ){
		vector<string> tokvec = tok.Tokenize(idvec[i], "\t");
		tid_g_map.insert(make_pair( tokvec[0],tokvec[1] ) ); 
	}
	
	vector<string> exprvec  = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/expr/binary/1.rc");
	for(int i =0; i<exprvec.size(); i++ ){
		vector<string> tokvec = tok.Tokenize(exprvec[i], "\t");
		exprmap.insert(make_pair( tokvec[0],atof(tokvec[1].c_str())) ); 
	}
	
	map<string, double> spmap;
	vector<string> spvec  = tio.dlread_vector(spfname);
	for(int i =0; i<spvec.size(); i++ ){
		
		spmap.insert(make_pair( spvec[0],0) ); 
	}
	
	
	vector<string> glistvec  = tio.dlread_vector(glistfname);
	for(int i =0; i<glistvec.size(); i++ ){
		
		vector<string> tokvec = tok.Tokenize( g_tid_map.find(glistvec[i])->second, "_");
		
		//for every transcript 
		for(int j =1 ;j<tokvec.size() ; j++){
			string tid = tokvec[j];
			
			ss.str("");
			double r = 1*((double) rand() / (RAND_MAX))+4;
			double rate = pow(2,r);
			
			ss << tid<<"\t" << exprmap.find(tid)->second/rate ;
			
			
			outVec.push_back(ss.str());
		}
		
	}
	
	
	return  outVec;
}

std::vector<std::string> Sim::binary_ds_RC(std::string glistfname,std::string spfname ){
	vector<string> outVec;
	
	TextIO tio;
	TOK tok;
	stringstream ss;
	ss.precision(15);
	map<string , string> tid_g_map ;
	map<string , string> g_tid_map =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );
	map<string, double> exprmap;
	vector<string> idvec  = tio.dlread_vector("/rhome/ywyang/bigdata/RefSeq/hg38/Transcript_Gene_ID");
	for(int i =0; i<idvec.size(); i++ ){
		vector<string> tokvec = tok.Tokenize(idvec[i], "\t");
		tid_g_map.insert(make_pair( tokvec[0],tokvec[1] ) ); 
	}
	
	vector<string> exprvec  = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/expr/binary/1.rc");
	for(int i =0; i<exprvec.size(); i++ ){
		vector<string> tokvec = tok.Tokenize(exprvec[i], "\t");
		exprmap.insert(make_pair( tokvec[0],atof(tokvec[1].c_str())) ); 
	}
	
	map<string, double> spmap;
	vector<string> spvec  = tio.dlread_vector(spfname);
	for(int i =0; i<spvec.size(); i++ ){
		
		spmap.insert(make_pair( spvec[0],0) ); 
	}
	
	
	vector<string> glistvec  = tio.dlread_vector(glistfname);
	for(int i =0; i<glistvec.size(); i++ ){
		
		vector<string> tokvec = tok.Tokenize( g_tid_map.find(glistvec[i])->second, "_");
		
		//for every transcript 
		double min = DBL_MAX;
			double max = -1;
			int min_idx =0 ;
		int max_idx =0 ;
		vector<double> rcvec ;
		for(int j =1 ;j<tokvec.size() ; j++){
			
			string tid = tokvec[j];
			double expr = exprmap.find(tid)->second;
			rcvec.push_back(expr) ;
			
			if(expr < min){
				min = expr;
				min_idx = j;
			}
			if(expr > max){
				max = expr;
				max_idx = j;
			}
		}
	
		//next_permutation(tokvec.begin(),tokvec.end());
		iter_swap(tokvec.begin() + min_idx, tokvec.begin() + max_idx);
		for(int j =1 ;j<tokvec.size() ; j++){
			ss.str("");
			ss<< tokvec[j]<< "\t" << rcvec[j-1];
			outVec.push_back(ss.str()) ;
		}
		
		
	}
	
	
	return  outVec;
}

std::vector<std::string> Sim::binary_ALL_RC(std::string glistfname, std::vector<int> GroupVec, double dispersion){
	TextIO tio;
	TOK tok;
	stringstream ss,ssout;
	vector<string> outVec;
	vector< vector<string> > rc1vec = tio.dlread("/rhome/ywyang/bigdata/AS/sim/expr/binary/1.rc", "\t");
	map<string, double> rc1map ;  
	for(int i =0 ; i<rc1vec.size() ; i++){
		rc1map.insert( make_pair( rc1vec[i][0],atof( rc1vec[i][1].c_str()) ) );
	}
	
	
	vector< vector<string> > rc2vec = tio.dlread("/rhome/ywyang/bigdata/AS/sim/expr/binary/2.rc", "\t");
	map<string, double> rc2map ;  
	for(int i =0 ; i<rc2vec.size() ; i++){
		rc2map.insert( make_pair( rc2vec[i][0], atof( rc2vec[i][1].c_str()) ) );
	}
	
	vector<map< string, double> > rcmapvec;
	rcmapvec.push_back(rc1map);
	rcmapvec.push_back(rc2map);
	
	
	
	map<string , string> g_tid_map =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );

	vector<string> glistvec = tio.dlread_vector(glistfname); 
	 
	for(int i =0  ; i<glistvec.size() ; i++){
		string gid  = glistvec[i];
		string tidstr = g_tid_map.find(gid)->second;
		vector<string> tidvec = tok.Tokenize( tidstr, "_" );
		//cout <<tidstr << endl;
			//for every transcript 
		for(int j =1 ; j< tidvec.size() ; j++ ){
			string tid = tidvec[j] ; 
			ssout.str("");
			//cout << GroupVec.size() << endl;
			for(int n = 0 ; n< GroupVec.size(); n++){
				map<string, double> & rcmap = rcmapvec[n];
				int sample_size = GroupVec[n];
				 
				//sample #sample_size read counts.
				ss.str("");
				double mean = rcmap.find(tid)->second; 
				ss << "Rscript RCGen.r "<< mean << " " << dispersion<< " " <<sample_size ;
				//cout << ss.str()<< endl; 
				string cmd = ss.str();
				system(cmd.c_str());
				
				//
				
				if(n==0){
					ssout << tid;
				}
				vector<string> tmpvec = tio.dlread_vector("rctmp");
 				for(int k= 0 ; k<sample_size; k++){
 					ssout << "\t"<< tmpvec[k];
 				}
 				
 				
			}
			
			outVec.push_back(ssout.str());
			
			
		}
		
	}
	
	
	return outVec;
}

std::vector<std::string> Sim::binary_ALL_RC(std::string glistfname, std::vector<int> GroupVec, double dispersion,int IDX){
	TextIO tio;
	TOK tok;
	stringstream ss,ssout;
	vector<string> outVec;
	vector< vector<string> > rc1vec = tio.dlread("/rhome/ywyang/bigdata/AS/sim/expr/binary/1.rc", "\t");
	map<string, double> rc1map ;  
	for(int i =0 ; i<rc1vec.size() ; i++){
		rc1map.insert( make_pair( rc1vec[i][0],atof( rc1vec[i][1].c_str()) ) );
	}
	
	
	vector< vector<string> > rc2vec = tio.dlread("/rhome/ywyang/bigdata/AS/sim/expr/binary/2.rc", "\t");
	map<string, double> rc2map ;  
	for(int i =0 ; i<rc2vec.size() ; i++){
		rc2map.insert( make_pair( rc2vec[i][0], atof( rc2vec[i][1].c_str()) ) );
	}
	
	vector<map< string, double> > rcmapvec;
	rcmapvec.push_back(rc1map);
	rcmapvec.push_back(rc2map);
	
	
	
	map<string , string> g_tid_map =  read_transcripts_map("/rhome/ywyang/bigdata/AS/sim/gene_list/all_genes.txt", "/rhome/ywyang/bigdata/RefSeq/hg38/gene_gtf/" );

	vector<string> glistvec = tio.dlread_vector(glistfname); 
	 
	for(int i =0  ; i<glistvec.size() ; i++){
		if(i%30 == IDX){	
	  cout << i << endl;
		string gid  = glistvec[i];
		string tidstr = g_tid_map.find(gid)->second;
		vector<string> tidvec = tok.Tokenize( tidstr, "_" );
		//cout <<tidstr << endl;
			//for every transcript 
		
			for(int j =1 ; j< tidvec.size() ; j++ ){
			string tid = tidvec[j] ; 
			ssout.str("");
			//cout << GroupVec.size() << endl;
			for(int n = 0 ; n< GroupVec.size(); n++){
				map<string, double> & rcmap = rcmapvec[n];
				int sample_size = GroupVec[n];
				 
				//sample #sample_size read counts.
				ss.str("");
				double mean = rcmap.find(tid)->second; 
				ss << "Rscript RCGen_para.r "<< mean << " " << dispersion<< " " <<sample_size<< " "<<IDX ;
				//cout << ss.str()<< endl; 
				string cmd = ss.str();
				system(cmd.c_str());
				
				//
				
				if(n==0){
					ssout << tid;
				}
				ss.str("");
				ss<< "rctmp"<<IDX;
				vector<string> tmpvec = tio.dlread_vector(ss.str());
 				for(int k= 0 ; k<sample_size; k++){
 					ssout << "\t"<< tmpvec[k];
 				}
 				
 				
			}
			
			outVec.push_back(ssout.str());
			
			
			}
		}
		
	}
	
	
	return outVec;
}


void Sim::naive_1(int size_a, int size_b){
	stringstream ss;
	TextIO tio;
	TOK tok;
	map<string,string> ExprMap,BedMap;
	//reads the expr file
	vector<string> exprvec = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/naive/gtf_beds/3eg.expr");
	for(int i =0; i<exprvec.size(); i++){
		vector<string> tmpvec = tok.Tokenize(exprvec[i],"\t");
		ExprMap.insert(make_pair(tmpvec[0],exprvec[i]));
	}
	//reads the bedfile
	vector<string> bedvec = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/naive/gtf_beds/3eg.bed");
	for(int i =0; i<bedvec.size(); i++){
		vector<string> tmpvec = tok.Tokenize(bedvec[i],"\t");
		BedMap.insert(make_pair(tmpvec[3],bedvec[i]));
	}
	
	
	
	
	for(int i =0  ; i<100 ; i++){
		vector<int> RCVec_t1,RCVec_t2;  
		//generate rc 
		//t1
		ss.str("");
		ss<< "Rscript RCGen.r " << 600 <<" " << 0.179 << " " << size_a ; 
		string cmd = ss.str();
		system(cmd.c_str());
		
		vector<string> rctmpvec = tio.dlread_vector("rctmp");
		for(int j =0; j<rctmpvec.size(); j++){
			RCVec_t1.push_back(atoi(rctmpvec[j].c_str()));
		}
		
		ss.str("");
		ss<< "Rscript RCGen.r " << 600 <<" " << 0.179 << " " << size_b ; 
		cmd = ss.str();
		system(cmd.c_str());
		
		rctmpvec = tio.dlread_vector("rctmp");
		for(int j =0; j<rctmpvec.size(); j++){
			RCVec_t1.push_back(atoi(rctmpvec[j].c_str()));
		}
		
		//t2
		ss.str("");
		ss<< "Rscript RCGen.r " <<  400<<" " << 0.179 << " " << size_a ; 
		cmd = ss.str();
		system(cmd.c_str());
		
		rctmpvec = tio.dlread_vector("rctmp");
		for(int j =0; j<rctmpvec.size(); j++){
			RCVec_t2.push_back(atoi(rctmpvec[j].c_str()));
		}
		
		ss.str("");
		ss<< "Rscript RCGen.r " << 400<<" " << 0.179 << " " << size_b ; 
		cmd = ss.str();
		system(cmd.c_str());
		
		rctmpvec = tio.dlread_vector("rctmp");
		for(int j =0; j<rctmpvec.size(); j++){
			RCVec_t2.push_back(atoi(rctmpvec[j].c_str()));
		}
		
		//create bed
	
		for(int j =0 ; j<RCVec_t2.size(); j++){
			vector<string> readvec;
			vector<string> tmp_reads;
			
			
			//sim t1
			tmp_reads=simReads( ExprMap["t1"], BedMap["t1"], "t1", RCVec_t1[j] , 50 ,0) ;
			//read t1
			for(int k=0; k<tmp_reads.size(); k++){
				readvec.push_back(tmp_reads[k]);
			}
			//sim t2
			tmp_reads=simReads( ExprMap["t2"], BedMap["t2"], "t2", RCVec_t2[j] , 50 ,0) ;
			//read t2
			for(int k=0; k<tmp_reads.size(); k++){
				readvec.push_back(tmp_reads[k]);
			}
			//write reads
			ss.str("");
			ss<< "/rhome/ywyang/bigdata/AS/sim/naive/"<<size_a<<"_"<<size_b<<"/3eg_1_" << i<< "_"<<j<<".bed"   ;
			tio.dlwrite( ss.str(), readvec);
		}
		
		
				
	}
	
	
	
}

void Sim::naive_2(int size_a, int size_b){
	stringstream ss;
	TextIO tio;
	TOK tok;
	map<string,string> ExprMap,BedMap;
	//reads the expr file
	vector<string> exprvec = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/naive/gtf_beds/3eg.expr");
	for(int i =0; i<exprvec.size(); i++){
		vector<string> tmpvec = tok.Tokenize(exprvec[i],"\t");
		ExprMap.insert(make_pair(tmpvec[0],exprvec[i]));
	}
	//reads the bedfile
	vector<string> bedvec = tio.dlread_vector("/rhome/ywyang/bigdata/AS/sim/naive/gtf_beds/3eg.bed");
	for(int i =0; i<bedvec.size(); i++){
		vector<string> tmpvec = tok.Tokenize(bedvec[i],"\t");
		BedMap.insert(make_pair(tmpvec[3],bedvec[i]));
	}
	
	
	
	
	for(int i =0  ; i<20 ; i++){
		vector<int> RCVec_t1,RCVec_t2;  
		//generate rc 
		//t1
		ss.str("");
		ss<< "Rscript RCGen.r " << 900 <<" " << 0.179 << " " << size_a ; 
		string cmd = ss.str();
		system(cmd.c_str());
		
		vector<string> rctmpvec = tio.dlread_vector("rctmp");
		for(int j =0; j<rctmpvec.size(); j++){
			RCVec_t1.push_back(atoi(rctmpvec[j].c_str()));
		}
		
		ss.str("");
		ss<< "Rscript RCGen.r " << 300 <<" " << 0.179 << " " << size_b ; 
		cmd = ss.str();
		system(cmd.c_str());
		
		rctmpvec = tio.dlread_vector("rctmp");
		for(int j =0; j<rctmpvec.size(); j++){
			RCVec_t1.push_back(atoi(rctmpvec[j].c_str()));
		}
		
		//t2
		ss.str("");
		ss<< "Rscript RCGen.r " << 200 <<" " << 0.179 << " " << size_a ; 
		cmd = ss.str();
		system(cmd.c_str());
		
		rctmpvec = tio.dlread_vector("rctmp");
		for(int j =0; j<rctmpvec.size(); j++){
			RCVec_t2.push_back(atoi(rctmpvec[j].c_str()));
		}
		
		ss.str("");
		ss<< "Rscript RCGen.r " << 600 <<" " << 0.179 << " " << size_b ; 
		cmd = ss.str();
		system(cmd.c_str());
		
		rctmpvec = tio.dlread_vector("rctmp");
		for(int j =0; j<rctmpvec.size(); j++){
			RCVec_t2.push_back(atoi(rctmpvec[j].c_str()));
		}
		
		//create bed
	
		for(int j =0 ; j<RCVec_t2.size(); j++){
			vector<string> readvec;
			vector<string> tmp_reads;
			
			
			//sim t1
			tmp_reads=simReads( ExprMap["t1"], BedMap["t1"], "t1", RCVec_t1[j] , 100 ,0) ;
			//read t1
			for(int k=0; k<tmp_reads.size(); k++){
				readvec.push_back(tmp_reads[k]);
			}
			//sim t2
			tmp_reads=simReads( ExprMap["t2"], BedMap["t2"], "t2", RCVec_t2[j] , 100 ,0) ;
			//read t2
			for(int k=0; k<tmp_reads.size(); k++){
				readvec.push_back(tmp_reads[k]);
			}
			//write reads
			ss.str("");
			ss<< "/rhome/ywyang/bigdata/AS/sim/naive/"<<size_a<<"_"<<size_b<<"/3eg_2_" << i<< "_"<<j<<".bed"   ;
			tio.dlwrite( ss.str(), readvec);
		}
		
		
				
	}
	
	
	
}