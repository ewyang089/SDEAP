#include "SDEAP.h"

using namespace std;




void SDEAP::cat_all(std::string outdir ){
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
			PredVec.push_back(tmpvec[j]);
		}
	}
	tio.dlwrite( outdir+"/all_pred" , PredVec);
	
}



vector<string> SDEAP::var(std::string glistfname, int sample_size, double min_mean, double max_cv2 , double min_cv2_rate){
	TextIO tio;
	TOK tok;
	vector<string> glistvec = tio.dlread_vector(glistfname);
	vector<string> PredVec;
	vector<string> TmpResVec;
	map<string, string> PredMap;
	map<string, double> ScoreMap;
	stringstream ss;
	
	for(int i=0 ; i<glistvec.size() ; i++){
		double mean =0;
		//cout << "mean:"<< glistvec[i] << endl;
		vector<string> tokvec = tok.Tokenize(glistvec[i],"\t");
		ScoreMap.insert(make_pair( tokvec[0]+tokvec[1], atof(tokvec[tokvec.size()-1].c_str())));		
		for(int j =0 ; j<sample_size; j++){
			 mean += atof(tokvec[j+2].c_str());
		}
		mean = mean/ (double) sample_size;
		double variance =0 ;
		
		for(int j =0 ; j<sample_size; j++){
				//cout << "tokvec[j+1]:"<< tokvec[j+2]<< endl;		
			  //cout << "@"<< atof(tokvec[j+2].c_str())-mean << endl;
			  variance += fabs(atof(tokvec[j+2].c_str())-mean)*fabs(atof(tokvec[j+2].c_str())-mean);
		}
		
		variance = variance / (double) sample_size;
		if(mean > min_mean ){
			double CV2 = variance/(mean*mean);
			ss.str("");
			ss <<tokvec[0] << "\t"<< tokvec[1]<<"\t" << mean<< "\t" <<variance << "\t"<<CV2;
			TmpResVec.push_back(ss.str());
		}
	}
	
	tio.dlwrite("tmpvar",TmpResVec);
	//run R
	
	ss.str("");
	ss<< "Rscript reg.r tmpvar" ;
	string cmd = ss.str();
	system(cmd.c_str());
		
	//refine the raw prediction
	
	vector<string> cv2predvec =  tio.dlread_vector("tmpvar_res"); 
	map<string, string> selectedIdMap;
	for(int i =0 ; i< cv2predvec.size(); i++){
		vector<string> tokvec = tok.Tokenize(cv2predvec[i],"\t");
		double predcv2= atof(tokvec[5].c_str());
		double cv2= atof(tokvec[4].c_str());
		double avg= atof(tokvec[2].c_str());
		
		if(  predcv2 < max_cv2 && cv2/predcv2> min_cv2_rate /*&& avg>4.5*/){
			
			selectedIdMap.insert(make_pair( tokvec[0]+tokvec[1],"" ));
			//cout << "select: " <<tokvec[0] << tokvec[1]<< endl;	
		}
		
	}
	
	vector<string> qpredvec =  tio.dlread_vector(glistfname); 
	
	for(int i =0 ; i< qpredvec.size(); i++){
		vector<string> tokvec = tok.Tokenize(qpredvec[i],"\t");
		//cout << "chk:"<< tokvec[0]<< tokvec[1]<< endl;
		if(selectedIdMap.find(tokvec[0]+tokvec[1])!=selectedIdMap.end()){
			PredVec.push_back(qpredvec[i]);
		}
		else{
			ss.str("");
			for(int j = 0 ; j<tokvec.size()-2; j++){
				ss<<tokvec[j]<<"\t";
			}
			
			if(ScoreMap[tokvec[0]+tokvec[1]]< 0.9999){
				ss<<0.8<<"\t"<<0.8;
			}
			else{
				ss<<0.9999<<"\t"<<0.9999;
			}			
			PredVec.push_back(ss.str());
			
			
		}
	}

	return PredVec;
}


std::vector<std::string> SDEAP::toGeneFeatures(std::string predfname){
	TextIO tio;
	TOK tok;
	vector<string> outvec ; 
	vector<string> featurevec = tio.dlread_vector(predfname);
	map<string,string> gstrmap;
	map<string,double> gidmap;
	stringstream ss;
	for(int i =0;  i< featurevec.size();i++){
		vector<string> tokvec = tok.Tokenize(featurevec[i],"\t");
		vector<string> tmpvec = tok.Tokenize(tokvec[0],".");
		double pvalue = atof(tokvec[tokvec.size()-1].c_str());
		string gname = tmpvec[0]; 
		map< string, double>::iterator mitr =  gidmap.find(gname);
		if(mitr!=gidmap.end() ){
			
			//update
			if(pvalue < mitr->second){
				mitr->second = pvalue;
				gstrmap[gname] =  featurevec[i];
			}
		}
		else{
			//insert a new gene
			gidmap.insert(make_pair(gname,pvalue));
			gstrmap.insert(make_pair(gname,featurevec[i]));
			
		}
	}
	
	for(map< string, double>::iterator outitr =gidmap.begin()   ; outitr != gidmap.end(); outitr++){
		if(outitr->second < 0.1){
			ss.str("");
			string predstr = gstrmap.find(outitr->first)->second;
			ss <<outitr->first << "\t"<< predstr;
			outvec.push_back(ss.str());
		}
	}
	
	return outvec;
}

std::vector<std::string> SDEAP::toPCAdata(std::string predfname){
	//normailze to mean 0 variance 1
	
	TextIO tio;
	TOK tok;
	vector<vector<string> > featurevec = tio.dlread(predfname,"\t");
  int sample_size =0;
	stringstream ss;
	vector<string> tmpvec; 
	for(int i =0 ; i< featurevec.size(); i++){
		ss.str("");
		sample_size = (featurevec[i].size()-5)/2;
		ss<< featurevec[i][0];
		for(int j =0 ; j<sample_size; j++){
			
			//ss<< "\t"<< log(atof(featurevec[i][j+3].c_str())+0.0001);
			ss<< "\t"<< featurevec[i][j+3];
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
		int pos_rc= 0;
		int neg_rc= 0;
		ss<< tokvec[0];
		for(int j =1;j<tokvec.size() ; j++){
			double x = atof(tokvec[j].c_str());
				ss<< "\t"<< (x-mean)/std;
//		
		}
			outvec.push_back(ss.str());
	
		
	} 
	
	return outvec; 
	
}




void SDEAP::sum_RCs(std::string glistfname,std::string outdir){
	TextIO tio;
	vector<string> glistvec = tio.dlread_vector(glistfname);
	vector<string> PredVec;
	
	for(int i=0 ; i<glistvec.size() ; i++){
		string gname = glistvec[i];
		vector<string> tmpvec = tio.dlread_vector(outdir+"/"+gname);
		for(int j=0 ; j< tmpvec.size(); j++){
			PredVec.push_back(gname+":"+tmpvec[j]);
		}
	}
	tio.dlwrite( outdir+"/all_pred" , PredVec);
	
}

void SDEAP::runSDEAP(std::vector<std::string> BedVec, std::string gname,std::string gtfdir  , std::string tmpdir , std::string outdir,std::vector< double> lib_sizevec , std::vector< int> lenvec){
		
	cout << "run_SDEAP ..." << gname << endl;
	TOK tok;
	stringstream ss;
	double scale =0; 
	vector<string> tmpdirvec;
	for(int n =0; n< BedVec.size() ; n++){
			vector<string> tokvec = tok.Tokenize( BedVec[n] , "/");
			string fname  = tokvec[tokvec.size()-1];
			
			ss.str("");
			ss<< tmpdir<<"/" <<  fname<<"_"<<n ;
			string cmd = "mkdir "  + ss.str() ;
			tmpdirvec.push_back(ss.str());
			//cout << cmd<<endl;
			system(cmd.c_str());
			cmd = "touch "  + BedVec[n];
			//cout << cmd<<endl;
			system(cmd.c_str());
	}
	
	
	
	vector<string> spgvec;
	vector<string> spg_rc_vec;

	vector<string> PredVec;
	

	for(int i =0 ; i< BedVec.size() ;  i++){
			double lib_size = lib_sizevec[i] ;
			scale = 0.5/lib_size;
			string bedfile = BedVec[i];
			cout << "reading ... "<< bedfile << endl;
			vector<string> tokvec = tok.Tokenize(bedfile,"/");
			string fname  = tokvec[tokvec.size()-1];
			//ss.str("");
			
			//ss<< tmpdir<<"/" <<   fname;
			string tmpfilepath = tmpdirvec[i];
			//tmpdirvec.push_back(tmpfilepath);
		
			//test ASMs with # of paths < 4
			gene g_spg;
			spgvec = g_spg.init_Hybrid( gtfdir + "/"+gname, bedfile , gname, tmpfilepath,lenvec[i],scale);
			//test ASMs with # of paths > 4 by using RCs
			gene g_spg_rc;
			spg_rc_vec = g_spg.init_Hybrid_RC(gtfdir +"/"+gname,  bedfile, gname, tmpfilepath,lenvec[i],scale);
	}

	
		for( int i = 0 ; i<spgvec.size() ; i++){
			//cout << "merge ...... " << spgvec[i] << endl;
			double lib_size = lib_sizevec[i] ;
			scale = 0.5/lib_size;
			Merge mg;
			//
			vector<string> resvec = mg.merge_Paths(tmpdirvec,outdir, spgvec[i],15*scale);	
			map<string, string> resmap;
			
			for(int j =0 ; j< resvec.size() ; j++){
			//cout << resvec[j]<< endl;
				PredVec.push_back(resvec[j]);
				vector<string> tokvec = tok.Tokenize(resvec[j],"\t");
				resmap.insert(make_pair(tokvec[1],resvec[j]));
			}
			
			//get gene features
			for(int j =0 ;j<mg.FNameVec.size(); j++){
			
				if(resmap.find(mg.FNameVec[j]) == resmap.end()){
					ss.str("");
					ss<< spgvec[i] << "\t";
					ss<< mg.FNameVec[j] << "\t";
					for(int k=0;k<mg.FeatureVec[j].size();k++){
						ss <<mg.FeatureVec[j][k] << "\t";
					}
					for(int k=0;k<mg.FeatureVec[j].size();k++){
						ss <<0 << "\t";
					}
					ss << 0.9999;
					resmap.insert(make_pair(mg.FNameVec[j],ss.str()));
					PredVec.push_back(ss.str());
				}
			}
			
			
		
			
		}

		for( int i = 0 ; i<spg_rc_vec.size() ; i++){
		  Merge mg;
		  double lib_size = lib_sizevec[i] ;
			scale = 0.5/lib_size;
			vector<string> resvec = mg.merge_RCs(tmpdirvec, outdir, spg_rc_vec[i],15*scale);	
			map<string, string> resmap;
			//cout << "merge ...... " << spg_rc_vec[i] << endl;
			for(int j =0 ; j< resvec.size() ; j++){
			//cout << resvec[j]<< endl;
				PredVec.push_back(resvec[j]);
				vector<string> tokvec = tok.Tokenize(resvec[j],"\t");
				resmap.insert(make_pair(tokvec[1],resvec[j]));
			}
			
			for(int j =0 ;j<mg.FNameVec.size(); j++){
				if(resmap.find(mg.FNameVec[j]) == resmap.end()){
					ss.str("");
					ss<< spgvec[i] << "\t";
					ss<< mg.FNameVec[j] << "\t";
					for(int k=0;k<mg.RCVec[j].size();k++){
						ss <<mg.RCVec[j][k] << "\t";
					}
					for(int k=0;k<mg.RCVec[j].size();k++){
						ss <<0 << "\t";
					}
					ss << 0.9999;
					resmap.insert(make_pair(mg.FNameVec[j],ss.str()));
					PredVec.push_back(ss.str());
				}
			}
			
		}
		
			
		TextIO tio;
		ss.str("");
		ss<<  outdir<<"/"<<gname;
		tio.dlwrite(ss.str(), PredVec);
		
		for(int n =0; n< tmpdirvec.size() ; n++){
			string cmd = "rm -rf "  + tmpdirvec[n] ;
			system(cmd.c_str());
		}
}

void SDEAP::runRC(std::vector<std::string> BedVec, std::string gname, std::string gtfdir , std::string tmpdir , std::string outdir,std::vector< double> lib_sizevec ){
		
	cout << "run_RC ..." << gname << endl;
	TOK tok;
	stringstream ss;
	double scale =0; 
	TextIO tio;
	vector<string> tmpdirvec;
	for(int n =0; n< BedVec.size() ; n++){
			vector<string> tokvec = tok.Tokenize( BedVec[n] , "/");
			string fname  = tokvec[tokvec.size()-1];
			
			ss.str("");
			ss<< tmpdir<<"/" <<  fname<<"_"<<n ;
			string cmd = "mkdir "  + ss.str() ;
			//cout << cmd<<endl;
			tmpdirvec.push_back(ss.str());
			system(cmd.c_str());
			cmd = "touch "  + BedVec[n];
			//cout << cmd<<endl;
			system(cmd.c_str());
	}
	
	
	
	vector<string> spgvec;

	vector<string> PredVec;
	

	for(int i =0 ; i< BedVec.size() ;  i++){
			double lib_size = lib_sizevec[i] ;
			scale = 0.5/lib_size;
			string bedfile = BedVec[i];
			cout << "reading ... "<< bedfile << endl;
			vector<string> tokvec = tok.Tokenize(bedfile,"/");
			string fname  = tokvec[tokvec.size()-1];
			
			string tmpfilepath = tmpdirvec[i];
			gene g_spg;
			spgvec = g_spg.init(gtfdir +"/"+gname, bedfile , gname, scale);
			ss.str("");
			ss<<  tmpfilepath<<"/"<<gname<<".rc";
			
			tio.dlwrite(ss.str(), spgvec);
	
	}

	Merge mg;
	
	PredVec = mg.mergeOneGene4Dexus( tmpdirvec,  gname);
			
		
	ss.str("");
	ss<<  outdir<<"/"<<gname;
	tio.dlwrite(ss.str(), PredVec);
		
		for(int n =0; n< BedVec.size() ; n++){
			vector<string> tokvec = tok.Tokenize( BedVec[n] , "/");
			string fname  = tokvec[tokvec.size()-1];
			
			string cmd = "rm -rf "  + tmpdirvec[n] ;
			system(cmd.c_str());
		}
}



int SDEAP::validCOV(std::string covfname){
	int outflag =0;
	int valid =0;
	TextIO tio;
	
	vector<vector<double> > depvec = tio.dlread_double(covfname, "\t");
	
	for(int i=0 ; i<depvec.size(); i++){
		int tmp_valid =1;
		vector<double> dvec = depvec[i];
		vector<double> covvec;
		//chk zeros propotions
		int zero_num =0;
		
		for(int j =0 ;j<dvec.size() ; j++){
			if(dvec[j]==0){
				 zero_num++;
			}
			else{
				covvec.push_back(dvec[j]);
			}
		}
		
		double zero_rate = zero_num/(double)dvec.size();
		
		if(zero_rate > 0.9){
			tmp_valid = 0;
		}
	
		//find the median of covered
		if(covvec.size()>0){
			sort(covvec.begin(), covvec.end());
			int median_idx = covvec.size()/2;
			if(covvec[median_idx] < 5){
				tmp_valid = 0;
			}
			

		}
		valid+=tmp_valid;
	}
	
	if(valid >=3){
		outflag =1;
	}	
	
	return outflag;
}
void SDEAP::runCOV(std::vector<std::string> BedVec, std::string gname,std::string gtfdir, std::string tmpdir , std::string outdir , std::vector< int> lenvec ){
		
	
	TOK tok;
	stringstream ss;
	
	GenomeSegment gs ; 
	gs.initGenomeSegment(gtfdir+"/"+gname, INT_MAX);
	//gs._debug();
	vector<string> covstrVec;
  cout << "run COVs...." << gname << endl; 
	for(int n =0; n< BedVec.size() ; n++){
			vector<string> tokvec = tok.Tokenize( BedVec[n] , "/");
			string fname  = tokvec[tokvec.size()-1];
			
			ss.str("");
			string cmd = "touch "  + BedVec[n];
			system(cmd.c_str());
	}
 
	for(int i =0;i< BedVec.size() ;i++ ){
	//
	vector<double> covvec =  gs.readscov(BedVec[i],lenvec[i]);
	
	ss.str("");
	for(int j =0 ;j<covvec.size() ; j++){
		ss << covvec[j] ;
		if(j < covvec.size()-1){
			ss<< "\t";
		}
	}	
	covstrVec.push_back(ss.str());
	//covstrVec.push_back(ss.str());
	
	
}

TextIO tio;
tio.dlwrite(outdir +"/"+gname+".cov",   covstrVec);


}

std::vector<double> SDEAP::calcSizeFactor(std::vector<std::string> fnamevec){
	TextIO tio;
	TOK tok;
	vector<double> sumvec;
	vector<double> outvec;
	for(int i =0 ; i<fnamevec.size(); i++){
		vector<string> sizevec = tio.dlread_vector(fnamevec[i]);
		double sum =0;
		for(int j =0 ; j<sizevec.size() ; j++){
			sum += atof(sizevec[j].c_str()); 
		}
		
		sumvec.push_back(sum);		
	}
	 
	for(int i =0 ; i<sumvec.size(); i++){
		outvec.push_back(sumvec[i]/1000000);

	}
	return outvec;
}

std::vector<double> SDEAP::Normalize4Dexus(std::vector<std::string> fnamevec){
	TextIO tio;
	TOK tok;
	vector<double> sumvec;
	vector<double> outvec;
	for(int i =0 ; i<fnamevec.size(); i++){
		
		vector<string> sizevec = tio.dlread_vector(fnamevec[i]);
		double sum =0;
		for(int j =0 ; j<sizevec.size() ; j++){
			sum += atof(sizevec[j].c_str()); 
		}
		
		sumvec.push_back(sum);		
	}
	double mean =0;
	for(int i =0 ; i<sumvec.size(); i++){
		mean+=sumvec[i];

	}
	mean = mean/sumvec.size(); 
	for(int i =0 ; i<sumvec.size(); i++){
		outvec.push_back(sumvec[i]/mean);

	}
	return outvec;
}

std::vector<std::string> 	SDEAP::toGeneRC(std::string indir, std::vector<std::string> gnamevec){
	vector<string> outVec;
	TextIO tio;
	TOK tok;
	stringstream ss;
	for(int i =0 ; i<gnamevec.size(); i++){
		vector< vector<string> > tokvec =  tio.dlread(indir+"/"+gnamevec[i],"\t"); 
		
		if(tokvec.size()>0){
			vector<double> dvec ;
			for(int j =0; j<tokvec.size(); j++ ){
				for(int k =1 ; k< tokvec[j].size(); k++){
					if(j==0){
						dvec.push_back(atof(tokvec[j][k].c_str() ));
					}
					else{
						dvec[k-1]+=atof(tokvec[j][k].c_str() );
					}
				}
			}
			ss.str("");
			ss<< gnamevec[i];
			for(int j=0; j<dvec.size() ; j++){
				ss<< "\t"<< dvec[j];
			}
			outVec.push_back(ss.str());
		}
		
	}
	
	return outVec;
}

