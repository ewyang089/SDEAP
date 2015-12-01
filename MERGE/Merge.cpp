#include "Merge.h"
using namespace std;

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

void Merge::merge4DexusExon(std::vector<std::string> pathvec, std::string outfname , std::string gene_list){
	vector<string> outVec;
	
	TextIO tio;
	stringstream ss;
	TOK tok;
	vector < vector<string> > gnamevec = tio.dlread( gene_list,"\t");
	
	for(int i =0; i<pathvec.size(); i++){
		ss<< "\t"<<i+1 ;
	}
	outVec.push_back(ss.str());
	
	for(int i =0; i<gnamevec.size(); i++){
		cout << gnamevec[i][0]<<endl;
		vector<string> tmpVec = mergeOneGene4Dexus(pathvec , gnamevec[i][0]);
		for(int j=0; j < tmpVec.size(); j++){
			outVec.push_back(gnamevec[i][0]+":"+tmpVec[j]);
		}	
	}
	
	tio.dlwrite(outfname, outVec);
	
}


void Merge::merge4DexusGene(std::vector<std::string> pathvec, std::string outfname , std::string gene_list){
	vector<string> outVec;
	
	TextIO tio;
	stringstream ss;
	TOK tok;
	vector < vector<string> > gnamevec = tio.dlread( gene_list,"\t");
	vector<double> grcvec ;
	for(int i =0; i<pathvec.size(); i++){
		ss<< "\t"<<i+1 ;
		grcvec.push_back(0.0);
	}
	outVec.push_back(ss.str());
	 
	for(int i =0; i<gnamevec.size(); i++){
		cout << gnamevec[i][0]<<endl;
		vector<string> tmpVec = mergeOneGene4Dexus(pathvec , gnamevec[i][0]);
		for(int j =0; j<grcvec.size(); j++){
					grcvec[j] = 0;
		}
		for(int j=0; j < tmpVec.size(); j++){
			
			vector<string> rcvec = tok.Tokenize(tmpVec[j], "\t");
			for(int k=1; k<rcvec.size(); k++){
					grcvec[k-1] += atof(rcvec[k].c_str());
			}
		}	
		ss.str("");
		ss<< gnamevec[i][0];
		for(int j =0; j<grcvec.size(); j++){
					ss<< "\t" << grcvec[j] ;
		}
		outVec.push_back(ss.str());
	}
	
	tio.dlwrite(outfname, outVec);
	
}


void Merge::merge_SPG(std::vector<std::string> pathvec, std::string outdir , std::string gene_list){
	
	
	TextIO tio;
	stringstream ss;
	TOK tok;
	vector < vector<string> > gnamevec = tio.dlread( gene_list,"\t");
	vector<string> resvec;
	for(int i =0; i<gnamevec.size(); i++){
		vector<string> outVec;
		cout << gnamevec[i][0]<<endl;
		vector<string> tmpVec = mergeOneGene_select(pathvec , gnamevec[i][0],2.0);
		if(tmpVec.size() > 0 ){
		for(int j=0; j < tmpVec.size(); j++){
			outVec.push_back(tmpVec[j]);
		}	
		
		tio.dlwrite(outdir+"/"+gnamevec[i][0], outVec);
		if(outVec.size()>0){
		if(outVec.size()==1){
			string cmd = "python DPMM_vector.py " +outdir+"/"+gnamevec[i][0];
			system(cmd.c_str());
		}
		if(outVec.size()>1){
			string cmd = "python DPMM.py " + outdir+"/"+gnamevec[i][0];
			system(cmd.c_str());
		}
		
		// delete  unreliable clustering
		vector<vector<string> >datavec = tio.dlread( outdir+"/"+gnamevec[i][0],"\t");
		vector<vector<string> >clustervec =  tio.dlread( outdir+"/"+gnamevec[i][0]+".cluster","\t");
		vector<vector<string> >datavec_new;
		vector<vector<string> >clustervec_new;
		
		for(int q  =0;q< clustervec.size(); q++){
			
			map<string,int>  clustermap,datamap;
			map<string,double> meanmap;
			int flag =0;
			for(int j=0; j<clustervec[q].size(); j++){
					if(clustermap.find(clustervec[q][j]) == clustermap.end()){
						clustermap.insert(make_pair( clustervec[q][j] , 0));
						datamap.insert(make_pair( clustervec[q][j] , 0));
						meanmap.insert(make_pair( clustervec[q][j] , 0));
					}
					
					(clustermap.find(clustervec[q][j])->second)++;
					double d = atof(datavec[q][j+1].c_str());
					if(d>1.0){
						(datamap.find(clustervec[q][j])->second)++;
					}
					meanmap.find(clustervec[q][j])->second +=d;
			}
			
			for(map<string, int >::iterator d_itr = datamap.begin() ; d_itr!= datamap.end(); d_itr++ ){
				
				if(d_itr->second>0){
					if(d_itr->second ==clustermap.find(d_itr->first)->second ){
						flag =1;
					}
				}
			}
			
			
		  /*
			//chk min and max mean of clusters
			double max =DBL_MIN; 
			double min=DBL_MAX ;
			for(map<string, double >::iterator d_itr = meanmap.begin() ; d_itr!= meanmap.end(); d_itr++ ){
				
				
				if(d_itr->second>0){
					
					d_itr->second = d_itr->second/(double)(clustermap.find(d_itr->first)->second);
					
				}
				
				if(d_itr->second > max){
					max = d_itr->second;
				}
				if(d_itr->second < min){
					min = d_itr->second;
				}
				
			}
			
			if(((max+0.01)/(min+0.01))<3){
				flag =0;
			}
			*/
			//
			
			if(flag ==1){
				datavec_new.push_back(datavec[q]);
				clustervec_new.push_back(clustervec[q]);
			} 			
		}
		
		if(datavec_new.size()>0){
		//output new data and cluster
		tio.dlwrite(outdir+"/"+gnamevec[i][0]+".data_cluster", clustervec_new,"\t");
		tio.dlwrite(outdir+"/"+gnamevec[i][0]+".data", datavec_new,"\t");
		// calcuate p-values by ANOVA
					
			string cmd = "Rscript ANOVA.r " + outdir+"/"+gnamevec[i][0]+".data "+ outdir+"/"+gnamevec[i][0]+".data_cluster ";
			cout << cmd<< endl;
			system(cmd.c_str());
		
			//merge p_value
			vector< string > fnamevec = tio.dlread_vector( outdir+"/"+gnamevec[i][0]+".data");
			vector< string > pvec = tio.dlread_vector( outdir+"/"+gnamevec[i][0]+".data.pvalue");
			
			for(int j =0 ; j<pvec.size() ; j++){
				resvec.push_back(gnamevec[i][0]+":"+fnamevec[j]+"\t"+pvec[j]);
			}
			
			
		
		}
		}
		}
	}
	tio.dlwrite(outdir+"/all_pred.txt", resvec);
	
}

std::vector<std::string> Merge::merge_Paths(std::vector<std::string> pathvec, std::string outdir , std::string gene_name,double cut){
	
	
	TextIO tio;
	stringstream ss;
	TOK tok;
	int min_clustersize = 1;
	vector<string> resvec;
	//for(int i =0; i<gnamevec.size(); i++){
		vector<string> outVec;
		//cout << gene_name<<endl;
		
		vector<string> tmpVec = mergeOnePath_select(pathvec , gene_name ,cut);
		
		if(tmpVec.size() > 0 ){
			for(int j=0; j < tmpVec.size(); j++){
				outVec.push_back(tmpVec[j]);
			}	
		
			tio.dlwrite(outdir+"/"+gene_name, outVec);
			if(outVec.size()>0){
				if(outVec.size()==1){
					string cmd = "python DPMM_vector.py " +outdir+"/"+gene_name;
					system(cmd.c_str());
					
				}
				if(outVec.size()>1){
					string cmd = "python DPMM.py " + outdir+"/"+gene_name;
					system(cmd.c_str());
					
				}
				
				
		// delete  unreliable clustering
				vector<vector<string> >datavec = tio.dlread( outdir+"/"+gene_name,"\t");
				vector<vector<string> >clustervec =  tio.dlread( outdir+"/"+gene_name+".cluster","\t");
				
				vector<vector<string> >datavec_new;
				vector<vector<string> >clustervec_new;
		
		for(int q  =0;q< clustervec.size(); q++){
			
					map<string,int>  clustermap,datamap;
					map<string,double> meanmap;
					int flag =0;
					for(int j=0; j<clustervec[q].size(); j++){
						if(clustermap.find(clustervec[q][j]) == clustermap.end()){
								clustermap.insert(make_pair( clustervec[q][j] , 0));
								datamap.insert(make_pair( clustervec[q][j] , 0));
								meanmap.insert(make_pair( clustervec[q][j] , 0));
						}
					
					(clustermap.find(clustervec[q][j])->second)++;
					double d = atof(datavec[q][j+1].c_str());
					if(d>1.0){
						(datamap.find(clustervec[q][j])->second)++;
					}
					meanmap.find(clustervec[q][j])->second +=d;
			}
			
			//chk if there is a cluster of features > 1.0
			for(map<string, int >::iterator d_itr = datamap.begin() ; d_itr!= datamap.end(); d_itr++ ){
				//cout << d_itr->first << "\t"<< d_itr->second << "\t" << (0.7* (double) clustermap.find(d_itr->first)->second) << endl;
				if(d_itr->second>0){
					if(d_itr->second >= 0.7* (double) clustermap.find(d_itr->first)->second ){
						
						flag =1;
					}
				}
			}
			// chk if the max mean is 4 times greater than the min mean
			double max_mean =-1 ;
			double min_mean = DBL_MAX;
				
			for(map<string, double >::iterator mean_itr = meanmap.begin() ; mean_itr !=meanmap.end() ; mean_itr ++){
				
				int size = clustermap.find(mean_itr->first)->second;
				double mean = mean_itr->second/ (double) size; 
				if(mean > max_mean){
					max_mean = mean;
				}
				if(mean < min_mean){
					min_mean = mean;
				}
			}
			//cout<<max_mean << ","<< min_mean<<","<<max_mean/min_mean<< endl;
			if((max_mean/min_mean) < 3.8/*2.5*/ ){
				flag =0;
			}
			//
		  

			if(flag ==1){
				
				datavec_new.push_back(datavec[q]);
				clustervec_new.push_back(clustervec[q]);
			} 			
		}
		
				if(datavec_new.size()>0){
		//output new data and cluster
					tio.dlwrite(outdir+"/"+gene_name+".data_cluster", clustervec_new,"\t");
					tio.dlwrite(outdir+"/"+gene_name+".data", datavec_new,"\t");
		// calcuate p-values by ANOVA
					
					string cmd = "Rscript ANOVA.r " + outdir+"/"+gene_name+".data "+ outdir+"/"+gene_name+".data_cluster ";
					//cout << cmd<< endl;
					system(cmd.c_str());
					//cmd = "rm -rf " + outdir+"/"+gene_name+".data "+ outdir+"/"+gene_name+".data_cluster ";
					//system(cmd.c_str());
		
			//merge p_value
					vector< vector< string > > fnamevec = tio.dlread( outdir+"/"+gene_name+".data", "\t");
					vector< string > pvec = tio.dlread_vector( outdir+"/"+gene_name+".data.pvalue");
			
					for(int j =0 ; j<pvec.size() ; j++){
						resvec.push_back(gene_name+"\t"+fnamevec[j][0]+"\t"+pvec[j]);
					}
					
					
					
					
						
				}
			}

		}
		string cmd;
		cmd = "rm -rf " +outdir+"/"+gene_name+".cluster " + outdir+"/"+gene_name;
		system(cmd.c_str());
					
		cmd = "rm -rf " + outdir+"/"+gene_name+".data "+ outdir+"/"+gene_name+".data_cluster " +outdir+"/"+gene_name+".data.pvalue";
		system(cmd.c_str());
		
		return resvec;
	
}

std::vector<std::string> Merge::merge_RCs(std::vector<std::string> pathvec, std::string outdir , std::string gene_name,double cut){
	
	
	TextIO tio;
	stringstream ss;
	TOK tok;

	vector<string> resvec;
	//for(int i =0; i<gnamevec.size(); i++){
		vector<string> outVec;
		//cout << gene_name<<endl;
		vector<string> tmpVec = mergeOneRC_select(pathvec , gene_name ,cut);
		
		
//		for(int i =0 ;i<tmpVec.size() ; i++){
//			
//			cout << " A : "<< tmpVec[i] << endl;
//		}
		
		if(tmpVec.size() > 0 ){
			for(int j=0; j < tmpVec.size(); j++){
				outVec.push_back(tmpVec[j]);
			}	
		
			tio.dlwrite(outdir+"/"+gene_name, outVec);
			if(outVec.size()>0){
				if(outVec.size()==1){
					string cmd = "python DPMM_vector.py " +outdir+"/"+gene_name;
					system(cmd.c_str());
				}
				if(outVec.size()>1){
					string cmd = "python DPMM.py " + outdir+"/"+gene_name;
					system(cmd.c_str());
				}
		// delete  unreliable clustering
				vector<vector<string> >datavec = tio.dlread( outdir+"/"+gene_name,"\t");
				vector<vector<string> >clustervec =  tio.dlread( outdir+"/"+gene_name+".cluster","\t");
				vector<vector<string> >datavec_new;
				vector<vector<string> >clustervec_new;
		
				for(int q  =0;q< clustervec.size(); q++){
			
			map<string,int>  clustermap,datamap;
			map<string,double> meanmap;
			int flag =0;
			for(int j=0; j<clustervec[q].size(); j++){
					if(clustermap.find(clustervec[q][j]) == clustermap.end()){
						clustermap.insert(make_pair( clustervec[q][j] , 0));
						datamap.insert(make_pair( clustervec[q][j] , 0));
						meanmap.insert(make_pair( clustervec[q][j] , 0));
					}
					
					(clustermap.find(clustervec[q][j])->second)++;
					double d = atof(datavec[q][j+1].c_str());
					if(d>1.0){
						(datamap.find(clustervec[q][j])->second)++;
					}
					meanmap.find(clustervec[q][j])->second +=d;
			}
			
			for(map<string, int >::iterator d_itr = datamap.begin() ; d_itr!= datamap.end(); d_itr++ ){
				
				if(d_itr->second>0){
					if(d_itr->second >= 0.7*clustermap.find(d_itr->first)->second ){
						flag =1;
					}
				}
			}
			double max_mean =-1 ;
			double min_mean = DBL_MAX;
				
			for(map<string, double >::iterator mean_itr = meanmap.begin() ; mean_itr !=meanmap.end() ; mean_itr ++){
				
				int size = clustermap.find(mean_itr->first)->second;
				double mean = mean_itr->second/ (double) size; 
				if(mean > max_mean){
					max_mean = mean;
				}
				if(mean < min_mean){
					min_mean = mean;
				}
			}
			//cout<<max_mean << ","<< min_mean<<","<<max_mean/min_mean<< endl;
			if((max_mean/min_mean) < 3.8/*2.5*/ ){
				flag =0;
			}
			
			if(flag ==1){
				datavec_new.push_back(datavec[q]);
				clustervec_new.push_back(clustervec[q]);
			} 			
				}
		
				if(datavec_new.size()>0){
		//output new data and cluster
					tio.dlwrite(outdir+"/"+gene_name+".data_cluster", clustervec_new,"\t");
					tio.dlwrite(outdir+"/"+gene_name+".data", datavec_new,"\t");
		// calcuate p-values by ANOVA
					
					string cmd = "Rscript ANOVA.r " + outdir+"/"+gene_name+".data "+ outdir+"/"+gene_name+".data_cluster ";
					
					system(cmd.c_str());
		
			//merge p_value
					vector< vector< string > > fnamevec = tio.dlread( outdir+"/"+gene_name+".data", "\t");
					vector< string > pvec = tio.dlread_vector( outdir+"/"+gene_name+".data.pvalue");
			
					for(int j =0 ; j<pvec.size() ; j++){
						resvec.push_back(gene_name+"\t"+fnamevec[j][0]+"\t"+pvec[j]);
					}
	
				}
			}

		}
		string cmd;
		cmd = "rm -rf " +outdir+"/"+gene_name+".cluster " + outdir+"/"+gene_name;
		system(cmd.c_str());
					
		cmd = "rm -rf " + outdir+"/"+gene_name+".data "+ outdir+"/"+gene_name+".data_cluster " +outdir+"/"+gene_name+".data.pvalue";
		system(cmd.c_str());
		return resvec;
	
}



std::vector<std::string> Merge::mergeOneGene(std::vector<std::string> pathvec , std::string gene_name){
	vector<string> outVec;
	vector<string> tmpVec;
	TextIO tio;
	stringstream ss;
	TOK tok;
	
	for(int i=0 ; i<pathvec.size(); i++){
		string fname = pathvec[i]+"/"+gene_name+".fpkm";
		if(!exists_test(fname.c_str())){
			cout << "Can't locate :"<<fname << endl;
			return outVec;
		}
			
		vector < vector<string> > tokvec = tio.dlread( fname,"\t");
		if(tokvec.size()==0){
			cout << fname <<" is empty." << endl;
			return outVec;
		}
	}
	
	
	for(int i=0 ; i<pathvec.size(); i++ ){
		string fname = pathvec[i]+"/"+gene_name+".fpkm";
		cout << "merge:" << fname <<endl;
		vector < vector<string> > tokvec = tio.dlread( fname,"\t");
		if(i==0){
			
			
			for(int j=0;j<tokvec.size();j++){
				ss.str("");
				ss<< tokvec[j][0] << "\t"<<tokvec[j][1];
				tmpVec.push_back(ss.str());
				
			}
			
		}
		else{
			for(int j=0;j<tokvec.size();j++){
				ss.str("");
				ss<<  "\t"<<tokvec[j][1];
				tmpVec[j]+=ss.str();
			}
		}
	}

	outVec = tmpVec;	
	return outVec;
}

std::vector<std::string> Merge::mergeOnePath_select(std::vector<std::string> pathvec , std::string gene_name,double cut){
	vector<string> outVec;
	vector<string> tmpVec;
	TextIO tio;
	stringstream ss;
	TOK tok;
	
	for(int i=0 ; i<pathvec.size(); i++){
		string fname = pathvec[i]+"/"+gene_name;
		
		if(!exists_test(fname.c_str())){
			cout << "Can't locate :"<<fname << endl;
			return outVec;
		}
			
		vector < vector<string> > tokvec = tio.dlread( fname,"\t");
		if(tokvec.size()==0){
			cout << fname <<" is empty." << endl;
			return outVec;
		}
	}

	vector<vector<int> > flagvec ;
	for(int i=0 ; i<pathvec.size(); i++){
		vector<vector<string> > tokvec = tio.dlread(pathvec[i]+"/"+gene_name,"\t");
		
		
		if(i ==0){
			for(int j =0 ;j<tokvec.size() ;j++){
				vector<int> tmpvec ;
				tmpvec.push_back(0); 
				flagvec.push_back(tmpvec);
			}
		}else{
			
			for(int j =0 ;j<tokvec.size() ;j++){
				
				flagvec[j].push_back(0);
				
			}
		}
		
	}
	
	
	
	for(int i=0 ; i<pathvec.size(); i++){
		vector<vector<string> > tokvec = tio.dlread(pathvec[i]+"/"+gene_name,"\t");
		
		
		if(i ==0){
			for(int j =0 ;j<tokvec.size() ;j++){
				vector<double> tmpvec ;
				tmpvec.push_back(atof(tokvec[j][1].c_str())); 
				FNameVec.push_back(tokvec[j][0]);
				FeatureVec.push_back(tmpvec);
			}
		}else{
			
			for(int j =0 ;j<tokvec.size() ;j++){
				double fpkm = atof(tokvec[j][1].c_str()); 
				FeatureVec[j].push_back(fpkm);
				
			}
		}
		
	}
	
	
	
	//calculate the RCs of paths
	
	for(int i=0 ; i<pathvec.size(); i++){
		vector<vector<string> > tokvec = tio.dlread(pathvec[i]+"/"+gene_name,"\t");
	
		if(i ==0){
			
			
			for(int j =0 ;j<tokvec.size() ;j++){
				//calc rc 
				double RC = 0;
				vector<string> nodevec = tok.Tokenize( tokvec[j][0], ",");
				for(int  k = 1 ; k < nodevec.size()-1;k++){
					vector<string> rcvec = tok.Tokenize( nodevec[k], ":");
					RC += atof(rcvec[2].c_str());
				}
				
				//
				
				vector<double> tmpvec ;
				tmpvec.push_back(RC); 
			
				RCVec.push_back(tmpvec);
			}
		}else{
			
			
			for(int j =0 ;j<tokvec.size() ;j++){
				//calc rc 
				double RC = 0;
				vector<string> nodevec = tok.Tokenize( tokvec[j][0], ",");
				for(int  k = 1 ; k < nodevec.size()-1;k++){
					vector<string> rcvec = tok.Tokenize( nodevec[k], ":");
					RC += atof(rcvec[2].c_str());
				}
				
				//
			
				RCVec[j].push_back(RC);
				
			}
		}
		
	}
	
	
	
	//
	/*
	for(int i =0; i < FeatureVec.size() ; i++){
		cout <<FNameVec[i] <<"\t";
		for(int j =0  ; j < FeatureVec[i].size(); j++){
			cout <<FeatureVec[i][j] <<"\t";
		}
		cout << endl;
	}
	cout << "RC  ====="<< endl;
	for(int i =0; i < FeatureVec.size() ; i++){
		cout <<FNameVec[i] <<"\t";
		for(int j =0  ; j < FeatureVec[i].size(); j++){
			cout <<RCVec[i][j] <<"\t";
		}
		cout << endl;
	}	
	*/	
	
// select major paths by %

	
	for(int i =0; i < FeatureVec.size() ; i++){
		vector<double> tmpvec;
		for(int j =0  ; j < FeatureVec[i].size(); j++){
			tmpvec.push_back(0);
		}
		percentvec.push_back(tmpvec);

	}	
	
	
	for(int j =0  ; j < FeatureVec[0].size(); j++){
		double col_sum =0;
		for(int i =0; i < FeatureVec.size() ; i++){
			col_sum += FeatureVec[i][j]; 
		}
		
		for(int i =0; i < FeatureVec.size() ; i++){
		
			percentvec[i][j] = FeatureVec[i][j]/col_sum; 
		}
		
	}	
//

	
	for(int i =0; i < FeatureVec.size() ; i++){
		int flag =0;
		for(int j =0  ; j < FeatureVec[i].size(); j++){
			if(RCVec[i][j]>cut && percentvec[i][j]>0.2 && FeatureVec[i][j] >1.0){
				flagvec[i][j] = 1;
			}
		}
		
	}
	
	for(int i =0; i < FeatureVec.size() ; i++){
		//cout <<FNameVec[i] <<"\t";
		ss.str("");
		int flag_sum = 0;
		for(int j =0  ; j < FeatureVec[i].size(); j++){
			flag_sum+=flagvec[i][j];
		}
		//cout << flag_sum  << endl;
		if(flag_sum>=2){
			ss<<FNameVec[i] << "\t";
			for(int j =0 ; j < FeatureVec[i].size(); j++){
				ss <<FeatureVec[i][j] ;
				if(j<FeatureVec[i].size()-1){
					ss << "\t";
				} 
			}
			outVec.push_back(ss.str());
		}
	
	}	
		
	return outVec;
}

std::vector<std::string> Merge::mergeOneRC_select(std::vector<std::string> pathvec , std::string gene_name,double cut){
	vector<string> outVec;
	vector<string> tmpVec;
	TextIO tio;
	stringstream ss;
	TOK tok;
	
	for(int i=0 ; i<pathvec.size(); i++){
		string fname = pathvec[i]+"/"+gene_name;
		if(!exists_test(fname.c_str())){
			cout << "Can't locate :"<<fname << endl;
			return outVec;
		}
			
		vector < vector<string> > tokvec = tio.dlread( fname,"\t");
		if(tokvec.size()==0){
			cout << fname <<" is empty." << endl;
			return outVec;
		}
	}
	
	vector< vector<double> > LenVec,RCVec, FeatureVec;
	vector< string > FNameVec;
	
  vector<vector<int> > flagvec ;
	
	for(int i=0 ; i<pathvec.size(); i++){
		vector<vector<string> > tokvec = tio.dlread(pathvec[i]+"/"+gene_name,"\t");
		
		if(i ==0){
			for(int j =0 ;j<tokvec.size() ;j++){
				vector<int> tmpvec ;
				tmpvec.push_back(0); 
				flagvec.push_back(tmpvec);
			}
		}else{
			
			for(int j =0 ;j<tokvec.size() ;j++){
				
				flagvec[j].push_back(0);
				
			}
		}
		
	}
	


	for(int i=0 ; i<pathvec.size(); i++){
		vector<vector<string> > tokvec = tio.dlread(pathvec[i]+"/"+gene_name,"\t");
	
		if(i ==0){
			for(int j =0 ;j<tokvec.size() ;j++){
				vector<double> tmpvec ;
				tmpvec.push_back(atof(tokvec[j][3].c_str())); 
				FNameVec.push_back(tokvec[j][0]);
				FeatureVec.push_back(tmpvec);
			}
		}
		else{
			
			for(int j =0 ;j<tokvec.size() ;j++){
				double fpkm = atof(tokvec[j][3].c_str()); 
				FeatureVec[j].push_back(fpkm);
				
			}
		}
		
	}
	
	


	for(int i=0 ; i<pathvec.size(); i++){
		vector<vector<string> > tokvec = tio.dlread(pathvec[i]+"/"+gene_name,"\t");
	
		if(i ==0){
			for(int j =0 ;j<tokvec.size() ;j++){
				vector<double> tmpvec ;
				tmpvec.push_back(atof(tokvec[j][2].c_str())); 
				RCVec.push_back(tmpvec);
			}
		}
		else{
			
			for(int j =0 ;j<tokvec.size() ;j++){
				double fpkm = atof(tokvec[j][2].c_str()); 
				RCVec[j].push_back(fpkm);
				
			}
		}
		
	}
	for(int i=0 ; i<pathvec.size(); i++){
		vector<vector<string> > tokvec = tio.dlread(pathvec[i]+"/"+gene_name,"\t");
	
		if(i ==0){
			for(int j =0 ;j<tokvec.size() ;j++){
				vector<double> tmpvec ;
				tmpvec.push_back(atof(tokvec[j][1].c_str())); 
				LenVec.push_back(tmpvec);
			}
		}
		else{
			
			for(int j =0 ;j<tokvec.size() ;j++){
				double fpkm = atof(tokvec[j][1].c_str()); 
				LenVec[j].push_back(fpkm);
				
			}
		}
		
	}

	
	for(int i =0; i < FeatureVec.size(); i++){
		for(int j =0; j < FeatureVec[0].size(); j++){
			
			//cout << flagvec[i][j]<< " "; 
			
			if(LenVec[i][j]>30 && RCVec[i][j] > cut && FeatureVec[i][j]> 1.0){
				flagvec[i][j]= 1;
			}
		}
		//cout << endl;
		
		int flag_sum =0; 
		for(int j =0; j < FeatureVec[0].size(); j++){
			flag_sum += flagvec[i][j];
		}
		ss.str("");
		if(flag_sum >=2){
			ss<<FNameVec[i] << "\t" ;
			for(int j =0; j < FeatureVec[0].size(); j++){
				ss<< FeatureVec[i][j];
				if(j < FeatureVec[0].size()-1){
					ss<< "\t";
				}
			}
			
			outVec.push_back(ss.str());
		}
	}
	
	
	return outVec;
}


std::vector<std::string> Merge::mergeOneGene_select(std::vector<std::string> pathvec , std::string gene_name,double cut){
	vector<string> outVec;
	vector<string> tmpVec;
	TextIO tio;
	stringstream ss;
	TOK tok;
	
	for(int i=0 ; i<pathvec.size(); i++){
		string fname = pathvec[i]+"/"+gene_name+".fpkm";
		if(!exists_test(fname.c_str())){
			cout << "Can't locate :"<<fname << endl;
			return outVec;
		}
			
		vector < vector<string> > tokvec = tio.dlread( fname,"\t");
		if(tokvec.size()==0){
			cout << fname <<" is empty." << endl;
			return outVec;
		}
	}
	
	map< string,string > JunctionMap;  
	for(int i=0 ; i<pathvec.size(); i++ ){
		string fname = pathvec[i]+"/"+gene_name+".fpkm";
	
		vector < vector<string> > tokvec = tio.dlread( fname,"\t");
						
		//if not junction....
		
		if(i==0){
			
			
			for(int j=0;j<tokvec.size();j++){
				vector<string> namevec = tok.Tokenize(tokvec[j][0],",");
				ss.str("");
   			ss<< tokvec[j][0] << "\t"<<tokvec[j][1];
				tmpVec.push_back(ss.str());
				
			
			}
			
		}
		else{
			for(int j=0;j<tokvec.size();j++){
				ss.str("");
				ss<<  "\t"<<tokvec[j][1];
				tmpVec[j]+=ss.str();
			}
		}
		
		
	}

	//calculate the mean of expression 
	vector<string> meanvec = tok.Tokenize(tmpVec[tmpVec.size()-1],"\t");
	double mean_sum = 0;
	for(int j =1 ; j<meanvec.size(); j++){
			mean_sum += atof(meanvec[j].c_str());
	}		
	
	double mean = mean_sum / (double)(meanvec.size()-1) ;
	//select non-zero rows 
	
	for(int i =0;i<tmpVec.size(); i++){
		
		vector<string> tokvec = tok.Tokenize(tmpVec[i],"\t");
		//calc the row mean
		double row_sum = 0;
		for(int j =1 ; j<tokvec.size(); j++){
			row_sum += atof(tokvec[j].c_str());
		}
		double row_mean = row_sum/(double)(tokvec.size()-1);
		
		int flag =0;
		for(int j =1 ; j<tokvec.size(); j++){
			double d = atof(tokvec[j].c_str());
			
			if(d>cut){
				flag =1;
			} 
		}
		if((row_mean/mean )<0.1){
			flag = 0;
		}
		if(flag ==1){
			//chk exon length
			vector<string> exonvec = tok.Tokenize(tokvec[0],",");
			
			if(exonvec.size() > 3){
				outVec.push_back(tmpVec[i]);
			}
			else{
				int exon_len=0;
				vector<string> nodevec = tok.Tokenize(exonvec[1],":");
				vector<string> lenvec = tok.Tokenize(nodevec[1],"$");
				exon_len = atoi(lenvec[1].c_str()) - atoi(lenvec[0].c_str()) +1;
				if(exon_len > 30){
					outVec.push_back(tmpVec[i]);
				}
			}
		}
		
	}
	
	return outVec;
}
std::vector<std::string> Merge::mergeOneGene4Dexus(std::vector<std::string> pathvec , std::string gene_name){
	vector<string> outVec;
	vector<string> tmpVec;
	TextIO tio;
	stringstream ss;
	TOK tok;
	
	for(int i=0 ; i<pathvec.size(); i++){
		string fname = pathvec[i]+"/"+gene_name+".rc";
		if(!exists_test(fname.c_str())){
			return outVec;
		}
			
		vector < vector<string> > tokvec = tio.dlread( fname,"\t");
		if(tokvec.size()==0){
			return outVec;
		}
	}
	
	
	for(int i=0 ; i<pathvec.size(); i++ ){
		string fname = pathvec[i]+"/"+gene_name+".rc";

		vector < vector<string> > tokvec = tio.dlread( fname,"\t");
		if(i==0){
			
			
			for(int j=0;j<tokvec.size();j++){
				ss.str("");
				ss<< tokvec[j][0] << "\t"<<tokvec[j][1];
				tmpVec.push_back(ss.str());
			}
			
		}
		else{
			for(int j=0;j<tokvec.size();j++){
				ss.str("");
				ss<<  "\t"<<tokvec[j][1];
				tmpVec[j]+=ss.str();
			}
		}
	}

	//select non-zero rows 
	for(int i =0;i<tmpVec.size(); i++){
		vector<string> tokvec = tok.Tokenize(tmpVec[i],"\t");
		int flag =0;
		//calculate fpkm
		vector<string> tmpvec = tok.Tokenize(tokvec[0],"$");
		double len = atof(tmpvec[1].c_str()) - atof(tmpvec[0].c_str()) +1;
		for(int j =1 ; j<tokvec.size(); j++){
			double d = atof(tokvec[j].c_str());
			double fpkm = d/len;
			if(fpkm>0){
				flag =1;
			} 
		}
		if(flag ==1){
			outVec.push_back(tmpVec[i]);
		}
		
	}
	return outVec;
}