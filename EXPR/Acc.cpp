#include "Acc.h"

using namespace std;
/*
void Acc::acc(){
	//ALL
	TextIO tio;
	map<string,double> AllMap;
	vector<vector<string> > allvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/ko_genes.txt", "\t"); 
	for(int i=0 ;i<allvec.size(); i++){
		AllMap.insert(make_pair(allvec[i][0] , 0.9999));	
	}
	//Real POS
	map<string,double> PosMap;
	vector<vector<string> > posvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/ko_pos.txt", "\t"); 
	for(int i=0 ;i<posvec.size(); i++){
		PosMap.insert(make_pair(posvec[i][0] , 0.9999));	
	}
	cout << PosMap.size()<< endl;
	//Real UP
	map<string,double> UpMap;
	vector<vector<string> > upvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/up_genes", "\t"); 
	for(int i=0 ;i<upvec.size(); i++){
		UpMap.insert(make_pair(upvec[i][0] , 0.9999));	
	}
	
	//Real Down
	map<string,double> DownMap;
	vector<vector<string> > downvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/down_genes", "\t"); 
	for(int i=0 ;i<downvec.size(); i++){
		DownMap.insert(make_pair(downvec[i][0] , 0.9999));	
	}
	
	//Real DS
	map<string,double> DSMap;
	vector<vector<string> > dsvec =  tio.dlread("/rhome/ywyang/bigdata/AS/sim/gene_list/ds_genes", "\t"); 
	for(int i=0 ;i<dsvec.size(); i++){
		DSMap.insert(make_pair(dsvec[i][0] , 0.9999));	
	}
	int TP,FP,TN,FN;
	TP=FP=TN=FN=0;
	//AS
	
	map<string,double> ASPredMap = readAS("/rhome/ywyang/bigdata/AS/sim/gene_list/ko_genes.txt", "/rhome/ywyang/bigdata/AS/sim/merge/SPG/7_3/all_pred.txt.qvalue");
	//cout << ASPredMap.size() << endl;
	
	double as_pre =0;
	double as_rec =0;
	
	

	for(map<string,double>::iterator as_itr  = ASPredMap.begin(); as_itr  != ASPredMap.end(); as_itr++){
		//cout << as_itr ->first << ","<< as_itr ->second<<endl;
		if(as_itr ->second <0.01){
			if(PosMap.find(as_itr ->first) != PosMap.end()){
				TP++;
			}
			else{
				FP++;
				cout << "FP:"<< as_itr ->first << "," << as_itr ->second<< endl;
			}
		}
		else{
			if(PosMap.find(as_itr ->first) != PosMap.end()){
				FN++;
				cout << "FN:"<< as_itr ->first << "," << as_itr ->second<< endl;
			}
			else{
				TN++;
			}
		}
	}
	
	cout << "AS ALL:"<< TP+FP<< ",PRE:" << TP/(double)(TP+FP)<< ",REC:" << TP/(double)(TP+FN)<< endl;

	//output 
	
  for(map<string,double>::iterator as_itr  = ASPredMap.begin(); as_itr  != ASPredMap.end(); as_itr++){
		cout << as_itr ->first << "\t"<< as_itr ->second<< "\t";
		
		cout<<endl;
	} 
	//exon Dexus
	
	map<string,double> DEPredMap = readDexus("/rhome/ywyang/bigdata/AS/sim/gene_list/ko_genes.txt", "/rhome/ywyang/bigdata/AS/sim/merge/Dexus/7_3/exon_pred.txt.dexus");
	//cout << ASPredMap.size() << endl;
	
	TP=FP=TN=FN=0;
	
	for(map<string,double>::iterator as_itr  = DEPredMap.begin(); as_itr  != DEPredMap.end(); as_itr++){
		//cout << as_itr ->first << ","<< as_itr ->second<<endl;
		if(as_itr ->second >0.5){
			if(PosMap.find(as_itr ->first) != PosMap.end()){
				TP++;
			}
			else{
				FP++;
				//cout << "FP:"<< as_itr ->first << "," << as_itr ->second<< endl;
			}
		}
		else{
			if(PosMap.find(as_itr ->first) != PosMap.end()){
				FN++;
				//cout << "FN:"<< as_itr ->first << "," << as_itr ->second<< endl;
			}
			else{
				TN++;
			}
		}
	}
	
	cout << "DE ALL:"<< TP+FP<< ",PRE:" << TP/(double)(TP+FP)<< ",REC:" << TP/(double)(TP+FN)<< endl;
	
}
*/


std::map<std::string, double> Acc::readAS(std::string allfname, std::string asfname){
	map<std::string, double>	AllMap;	
	map<std::string, string>	AllStrMap;	
	stringstream ss;
	TextIO tio;
	
	TOK tok;
	vector<vector<string> > allvec =  tio.dlread(allfname, "\t"); 
	vector<vector<string> > asvec =  tio.dlread(asfname, "\t");
	//cout <<allvec.size() << "," << asvec.size()<< endl;
	for(int i=0 ;i<allvec.size(); i++){
		AllMap.insert(make_pair(allvec[i][0] , 0.9999));
		AllStrMap.insert(make_pair(allvec[i][0] , ""));
		//cout << "add:"<< allvec[i][0]<<endl;	
	}	
	for(int i=0 ;i<asvec.size(); i++){
		vector<string> tokvec = tok.Tokenize(asvec[i][1],".");
		string gname = tokvec[0];
		map<string, double>::iterator all_itr = AllMap.find(gname) ;
		if(all_itr != AllMap.end()){
			if(all_itr->second > atof(asvec[i][asvec[i].size()-1].c_str() )){
				
				all_itr->second = atof(asvec[i][asvec[i].size()-1].c_str() );
				
				ss.str("");
				for(int j =0 ; j<asvec[i].size(); j++){
					ss<<asvec[i][j]; 
					if(j<asvec[i].size()-1){
						ss<<"\t";
					}
				}
				AllStrMap[gname]=ss.str(); 
				
			}	
		}
	}	
	
	return AllMap;
}

void Acc::readDexus(std::string out_dir ){

	map<std::string, double>	AllMap;	
	TextIO tio;
	TOK tok;
	stringstream ss;
	
	vector<vector<string> > predvec =  tio.dlread(out_dir +"/all_pred.dexus", "\t"); 
	vector<string > labelvec =  tio.dlread_vector(out_dir +"/all_pred.dexus_labels");
	//joint two files
	vector<string> mergedvec ;
	for(int i =0 ; i<predvec.size(); i++){
		ss.str("");
		
		for(int j =0; j<predvec[i].size()-2; j++){
			ss<<predvec[i][j] << "\t";
		}
		vector<string> tokvec = tok.Tokenize(labelvec[i] ,"\t" );
		for(int j =1 ; j<tokvec.size() ; j++){
			ss << tokvec[j]<< "\t";
		}
		
		ss<< predvec[i][predvec[i].size()-2]<< "\t"<< predvec[i][predvec[i].size()-1];
	  mergedvec.push_back(ss.str());
	}
	tio.dlwrite(out_dir+"/merged_pred.txt", mergedvec);
	//select gene features
	map<string, double> allmap;
	map<string, string> allstrmap;
	for(int i =0; i<mergedvec.size(); i++){
		vector<string> tokvec = tok.Tokenize( mergedvec[i],"\t");
		vector<string> tmpvec = tok.Tokenize( tokvec[0],":");
		string gname = tmpvec[0];
		double ini = atof(tokvec[tokvec.size()-2].c_str());
		map<string , double>::iterator mitr = allmap.find(gname);
		
		if(mitr == allmap.end() ){
			allmap.insert(make_pair( gname, ini));
			allstrmap.insert(make_pair( gname,mergedvec[i] ));
		}
		else{
			//update		
			if(ini  > mitr->second  ){
				allmap[gname] = ini;
				allstrmap[gname] = mergedvec[i];
			}
		}
	}
	vector<string> genevec;
	for(map<string,string>::iterator stritr = allstrmap.begin() ;stritr != allstrmap.end() ; stritr++ ){
		ss.str("");
		ss<<stritr->first << "\t"<<	stritr->second ;
		genevec.push_back(ss.str());	
	}
	
	tio.dlwrite(out_dir+"/gene_pred.txt", genevec);
	
	//sorting the gene features by the ini scores
	
	vector<double> inivec;
	for(int i =0 ; i< genevec.size(); i++){
		vector<string> tokvec = tok.Tokenize(genevec[i],"\t");
		double a = atof(tokvec[tokvec.size()-2].c_str());
		inivec.push_back(a);
	}
	
	
	for(int i =0 ; i< genevec.size(); i++){
		for(int j= 0; j <  genevec.size()-i-1  ; j++){
			
			double a = inivec[j];
			double b = inivec[j+1];
			if(a < b){
				//swap(a,b) in 
				double d= inivec[j];
				inivec[j] = inivec[j+1];
				inivec[j+1] = d;
				
				//swap str a and str b
				string tmpstr = genevec[j];
				genevec[j]= genevec[j+1];
				genevec[j+1] = tmpstr;
				 
			}
		}
	}
	
	 tio.dlwrite(out_dir+"/gene_sorted.txt", genevec);
	 
	 
	 
	
}
void Acc::readDexus_genes(std::string out_dir ){

	map<std::string, double>	AllMap;	
	TextIO tio;
	TOK tok;
	stringstream ss;
	
	vector<vector<string> > predvec =  tio.dlread(out_dir +"/gene.txt.dexus", "\t"); 
	vector<string > labelvec =  tio.dlread_vector(out_dir +"/gene.txt.dexus_labels");
	//joint two files
	vector<string> mergedvec ;
	for(int i =0 ; i<predvec.size(); i++){
		ss.str("");
		
		for(int j =0; j<predvec[i].size()-2; j++){
			ss<<predvec[i][j] << "\t";
		}
		vector<string> tokvec = tok.Tokenize(labelvec[i] ,"\t" );
		for(int j =1 ; j<tokvec.size() ; j++){
			ss << tokvec[j]<< "\t";
		}
		
		ss<< predvec[i][predvec[i].size()-2]<< "\t"<< predvec[i][predvec[i].size()-1];
	  mergedvec.push_back(ss.str());
	}
//	tio.dlwrite(out_dir+"/merged_pred.txt", mergedvec);
//	//select gene features
//	map<string, double> allmap;
//	map<string, string> allstrmap;
//	for(int i =0; i<mergedvec.size(); i++){
//		vector<string> tokvec = tok.Tokenize( mergedvec[i],"\t");
//		vector<string> tmpvec = tok.Tokenize( tokvec[0],":");
//		string gname = tmpvec[0];
//		double ini = atof(tokvec[tokvec.size()-2].c_str());
//		map<string , double>::iterator mitr = allmap.find(gname);
//		
//		if(mitr == allmap.end() ){
//			allmap.insert(make_pair( gname, ini));
//			allstrmap.insert(make_pair( gname,mergedvec[i] ));
//		}
//		else{
//			//update		
//			if(ini  > mitr->second  ){
//				allmap[gname] = ini;
//				allstrmap[gname] = mergedvec[i];
//			}
//		}
//	}
	vector<string> genevec = mergedvec;
//	for(map<string,string>::iterator stritr = allstrmap.begin() ;stritr != allstrmap.end() ; stritr++ ){
//		ss.str("");
//		ss<<stritr->first << "\t"<<	stritr->second ;
//		genevec.push_back(ss.str());	
//	}
//	
//	tio.dlwrite(out_dir+"/gene_pred.txt", genevec);
//	
	//sorting the gene features by the ini scores
	
	vector<double> inivec;
	for(int i =0 ; i< genevec.size(); i++){
		vector<string> tokvec = tok.Tokenize(genevec[i],"\t");
		double a = atof(tokvec[tokvec.size()-2].c_str());
		inivec.push_back(a);
	}
	
	
	for(int i =0 ; i< genevec.size(); i++){
		for(int j= 0; j <  genevec.size()-i-1  ; j++){
			
			double a = inivec[j];
			double b = inivec[j+1];
			if(a < b){
				//swap(a,b) in 
				double d= inivec[j];
				inivec[j] = inivec[j+1];
				inivec[j+1] = d;
				
				//swap str a and str b
				string tmpstr = genevec[j];
				genevec[j]= genevec[j+1];
				genevec[j+1] = tmpstr;
				 
			}
		}
	}
	
	 tio.dlwrite(out_dir+"/genelevel_sorted.txt", genevec);
	 
	 
	 
	
}

std::vector<std::string> Acc::toPCA(std::string out_dir, int g_num){
	
	TextIO tio;
	TOK tok;
	vector<vector<string> > InVec = tio.dlread(out_dir+"/gene_sorted.txt","\t");
	
	vector<vector<string> > featurevec;
	for(int i= 0 ; i<g_num; i++){
		featurevec.push_back(InVec[i]);
	}
	
	stringstream ss;
	vector<string> tmpvec; 
	for(int i =0 ; i< featurevec.size(); i++){
		ss.str("");
		int sample_size = (featurevec[i].size()-4)/2;
		ss<< featurevec[i][0];
		for(int j =0 ; j<sample_size; j++){
			ss<< "\t"<< featurevec[i][j+2];
		}
		tmpvec.push_back(ss.str());
	}
	
	vector<string> outvec;
	for(int i=0; i < tmpvec.size();i++){
		vector<string> tokvec = tok.Tokenize(tmpvec[i],"\t");
		//calc mean
		double mean=0;
		
		for(int j =1; j<tokvec.size(); j++ ){
			mean+= atof(tokvec[j].c_str());
		}
		mean = mean /(double)(tokvec.size()-1);
		//calc std
		double std=0;
		double sum = 0;
		for(int j =1; j<tokvec.size(); j++){
			double x = atof(tokvec[j].c_str());
			sum += (x-mean)*(x-mean);
			
		}
		sum = sum/(double)(tokvec.size()-1);
		std = sqrt(sum);
		//normalize
		ss.str("");
		
		ss<< tokvec[0];
		for(int j =1;j<tokvec.size() ; j++){
			double x = atof(tokvec[j].c_str());
			ss<< "\t"<< (x-mean)/std;
		}
		outvec.push_back(ss.str());
	} 
	
	return outvec; 
}

std::vector<std::string> Acc::toPCA_gene(std::string out_dir, int g_num){
	
	TextIO tio;
	TOK tok;
	vector<vector<string> > InVec = tio.dlread(out_dir+"/genelevel_sorted.txt","\t");
	
	vector<vector<string> > featurevec;
	for(int i= 0 ; i<g_num; i++){
		featurevec.push_back(InVec[i]);
	}
	
	stringstream ss;
	vector<string> tmpvec; 
	for(int i =0 ; i< featurevec.size(); i++){
		ss.str("");
		int sample_size = (featurevec[i].size()-3)/2;
		ss<< featurevec[i][0];
		for(int j =0 ; j<sample_size; j++){
			ss<< "\t"<< featurevec[i][j+1];
		}
		tmpvec.push_back(ss.str());
	}
	
	vector<string> outvec;
	for(int i=0; i < tmpvec.size();i++){
		vector<string> tokvec = tok.Tokenize(tmpvec[i],"\t");
		//calc mean
		double mean=0;
		
		for(int j =1; j<tokvec.size(); j++ ){
			mean+= atof(tokvec[j].c_str());
		}
		mean = mean /(double)(tokvec.size()-1);
		//calc std
		double std=0;
		double sum = 0;
		for(int j =1; j<tokvec.size(); j++){
			double x = atof(tokvec[j].c_str());
			sum += (x-mean)*(x-mean);
			
		}
		sum = sum/(double)(tokvec.size()-1);
		std = sqrt(sum);
		//normalize
		ss.str("");
		
		ss<< tokvec[0];
		for(int j =1;j<tokvec.size() ; j++){
			double x = atof(tokvec[j].c_str());
			ss<< "\t"<< (x-mean)/std;
		}
		outvec.push_back(ss.str());
	} 
	
	return outvec; 
}

//


std::vector<std::string>	Acc::readCuffdiff(std::string fname){
 	TextIO tio;
 	TOK tok ;
 	vector<string> outVec;
 	map<string,string> outMap;
 	vector<string> cuffvec = tio.dlread_vector(fname); 
 	for(int i =1;  i<cuffvec.size(); i++){
 		vector<string> tmpvec = tok.Tokenize(cuffvec[i], "\t");
 		/*
 		if(tmpvec[tmpvec.size()-1].compare("yes")==0){
 			if(outMap.find(tmpvec[1])==outMap.end()){
 				outVec.push_back(tmpvec[1]);
 				outMap.insert(make_pair(  tmpvec[1]  , "" ));
 			}
 		}*/
 		double fdr = atof(tmpvec[tmpvec.size()-2].c_str());
 		double lfr = atof(tmpvec[9].c_str());
 		lfr=fabs(lfr);
 		if(fdr<=0.2){
 			if(outMap.find(tmpvec[1])==outMap.end()){
 				outVec.push_back(tmpvec[1]);
 				outMap.insert(make_pair(  tmpvec[1]  , "" ));
 			}
 		}
 		
 	}
 	
 	return outVec;
}

void Acc::toDist(std::string corrfname){
	TextIO tio;
	vector<vector<double> > corrvec =  tio.dlread_double(corrfname,  "\t");
	vector<vector<double> > distmtx;
	
	for(int i =0 ; i<corrvec.size(); i++){
		vector<double> tmp_dist ; 
		for(int j =0 ; j<corrvec.size(); j++ ){
			double d = dist(corrvec[i][0],corrvec[i][1],corrvec[i][2], corrvec[j][0],corrvec[j][1],corrvec[j][2]);
			tmp_dist.push_back(d);
		}
		distmtx.push_back(tmp_dist);
	}
	tio.dlwrite_double( corrfname+".dist",distmtx,"\t");
}

double Acc::dist(double x1, double y1,double z1, double x2, double y2,double z2){
	
	return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)+ (z1-z2)*(z1-z2) );
}

void Acc::dexus_cat_all(std::string outdir ){
	TextIO tio;
	vector<string> glistvec; 
	
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (outdir.c_str())) != NULL) {
  /* print all the files and directories within directory */
  	while ((ent = readdir (dir)) != NULL) {
    	
    	string fname(ent->d_name);
    	if(fname.compare(".")!=0 && fname.compare("..")!=0){
    		glistvec.push_back(fname);
    	}
  	}
  	closedir (dir);
	} else {
  /* could not open directory */
  	perror ("");
  
	}
	
	vector<string> PredVec;
	
	for(int i=0 ; i<glistvec.size() ; i++){
		string gname = glistvec[i];
		if(i%100 == 0 ){
			cout << i << endl;
		}
		vector<string> tmpvec = tio.dlread_vector(outdir+"/"+gname);
		for(int j=0 ; j< tmpvec.size(); j++){
			PredVec.push_back(gname+":"+tmpvec[j]);
		}
	}
	tio.dlwrite( outdir+"/all_pred" , PredVec);
	
}


void Acc::to_pROC_sdeap(std::string predfname, std::string allfname, std::string dsfname,string outfname){
	TextIO tio;
	map<string,double> PredMap;
	map<string,int> LabelMap;
	stringstream ss;
	
	
	vector<string> allvec = tio.dlread_vector(allfname);
	
	//init maps
	
	for(int i =0; i<allvec.size();i++ ){
			PredMap.insert(make_pair(allvec[i] , 0.999));
			LabelMap.insert(make_pair(allvec[i], 0));
	}
	
	
	vector< vector<string> > predvec = tio.dlread(predfname , "\t");
	for(int i =0; i<predvec.size();i++){
		vector<string> ivec = predvec[i];
		string name = ivec[0];
		double score = atof(ivec[ivec.size()-1].c_str());
		if(PredMap.find(name)!= PredMap.end()){
			PredMap[name] =  score;
		}
	}
	
	vector<string> labelvec = tio.dlread_vector(dsfname);
	for(int i =0; i<labelvec.size();i++ ){
		string name =labelvec[i];
		if(LabelMap.find(name)!= LabelMap.end()){
			LabelMap[name] =  1;
		}
	}
	
	vector<string> outvec;
	for(int i =0; i<allvec.size();i++ ){
			ss.str("");
			string name = allvec[i];
			ss<< name<<"\t" << PredMap[name]<<"\t" << LabelMap[name];
			outvec.push_back(ss.str());
			
	}
	tio.dlwrite(outfname, outvec);
}

void Acc::to_pROC_dexus(std::string predfname, std::string allfname, std::string dsfname,std::string outfname){
	
	TextIO tio;
	map<string,double> PredMap;
	map<string,int> LabelMap;
	stringstream ss;
	
	
	vector<string> allvec = tio.dlread_vector(allfname);
	
	//init maps
	
	for(int i =0; i<allvec.size();i++ ){
			PredMap.insert(make_pair(allvec[i] , 0));
			LabelMap.insert(make_pair(allvec[i], 0));
	}
	
	
	vector< vector<string> > predvec = tio.dlread(predfname , "\t");
	for(int i =0; i<predvec.size();i++){
		vector<string> ivec = predvec[i];
		string name = ivec[0];
		double score = atof(ivec[ivec.size()-2].c_str());
		if(PredMap.find(name)!= PredMap.end()){
			PredMap[name] =  score;
		}
	}
	
	vector<string> labelvec = tio.dlread_vector(dsfname);
	for(int i =0; i<labelvec.size();i++ ){
		string name =labelvec[i];
		if(LabelMap.find(name)!= LabelMap.end()){
			LabelMap[name] =  1;
		}
	}
	
	vector<string> outvec;
	for(int i =0; i<allvec.size();i++ ){
			ss.str("");
			string name = allvec[i];
			ss<< name<<"\t" << PredMap[name]<<"\t" << LabelMap[name];
			outvec.push_back(ss.str());
			
	}
	tio.dlwrite(outfname, outvec);
}


void Acc::to_PR_sdeap(std::string predfname, std::string allfname, std::string dsfname,string outfname){
	TextIO tio;
	map<string,double> PredMap;
	map<string,int> LabelMap;
	stringstream ss;
	
	
	vector<string> allvec = tio.dlread_vector(allfname);
	
	//init maps
	
	for(int i =0; i<allvec.size();i++ ){
			PredMap.insert(make_pair(allvec[i] , 0.999));
			LabelMap.insert(make_pair(allvec[i], 0));
	}
	
	
	vector< vector<string> > predvec = tio.dlread(predfname , "\t");
	for(int i =0; i<predvec.size();i++){
		vector<string> ivec = predvec[i];
		string name = ivec[0];
		double score = atof(ivec[ivec.size()-1].c_str());
		if(PredMap.find(name)!= PredMap.end()){
			PredMap[name] =  score;
		}
	}
	
	vector<string> labelvec = tio.dlread_vector(dsfname);
	for(int i =0; i<labelvec.size();i++ ){
		string name =labelvec[i];
		if(LabelMap.find(name)!= LabelMap.end()){
			LabelMap[name] =  1;
		}
	}
	
	vector<string> outvec1,outvec2; 
	int TP,PredP;
	for(int i =0; i<allvec.size();i++ ){
			ss.str("");
			string name = allvec[i];
			ss<< name<<"\t" << PredMap[name]<<"\t" << LabelMap[name];
			
			
			if(LabelMap[name]==0){
			
				outvec1.push_back(ss.str());
				
				
			}
			else{
				
				outvec2.push_back(ss.str());
				
				
			}
	}
	tio.dlwrite(outfname+".0", outvec1);
	tio.dlwrite(outfname+".1", outvec2);
}


void Acc::to_PR_dexus(std::string predfname, std::string allfname, std::string dsfname,std::string outfname){
	
	TextIO tio;
	map<string,double> PredMap;
	map<string,int> LabelMap;
	stringstream ss;
	
	
	vector<string> allvec = tio.dlread_vector(allfname);
	
	//init maps
	
	for(int i =0; i<allvec.size();i++ ){
			PredMap.insert(make_pair(allvec[i] , 0));
			LabelMap.insert(make_pair(allvec[i], 0));
	}
	
	
	vector< vector<string> > predvec = tio.dlread(predfname , "\t");
	for(int i =0; i<predvec.size();i++){
		vector<string> ivec = predvec[i];
		string name = ivec[0];
		double score = atof(ivec[ivec.size()-2].c_str());
		if(PredMap.find(name)!= PredMap.end()){
			PredMap[name] =  score;
		}
	}
	
	vector<string> labelvec = tio.dlread_vector(dsfname);
	for(int i =0; i<labelvec.size();i++ ){
		string name =labelvec[i];
		if(LabelMap.find(name)!= LabelMap.end()){
			LabelMap[name] =  1;
		}
	}
	
	vector<string> outvec1,outvec2; 
	for(int i =0; i<allvec.size();i++ ){
			ss.str("");
			string name = allvec[i];
			ss<< name<<"\t" << PredMap[name]<<"\t" << LabelMap[name];
			if(LabelMap[name]==0){
				outvec1.push_back(ss.str());
			}
			else{
				outvec2.push_back(ss.str());
			}
	}
	tio.dlwrite(outfname+".0", outvec1);
	tio.dlwrite(outfname+".1", outvec2);
}