#include "Chrm.h"

using namespace std;

void Chrm::initChrm(string fname,int chr_len){
 gs.initGenomeSegment(fname, chr_len);
 
 //gs.printReads("../data/"); 
}




//void Chrm::createGeneClusters(std::string fname){
//	
//	vector<vector<int> > GClusterMtx;
//	//read gene number
//	int gnum = 0;
//	map<string, map<string, string> >::iterator exon_itr; 
//	map< string, int> GIdxmap;
//	 	
//	for(int i =0 ; i < gnum ; i++){
//		vector<int>
//		for(int j =0;j< gnum; j++){
//		}
//	}
//	
//	
//}


void Chrm::mapReads_BED(std::string fname){

	gs.mapReads_BED(fname);
}

int	Chrm::_check_split(std::string CICAR){
	int out_int = 0;
	
	for(int i=0;i<CICAR.length();i++ ){
		
		if(CICAR[i]=='N'){
			out_int =1;
			break;
		}
	}
	return out_int; 
}

void Chrm::printReadCounts(std::string outfname){
	gs.printReads(outfname); 
	
}
