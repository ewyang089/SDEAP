#include "Merge.h"
#include "TextIO.h"
#include "TOK.h"
#include "RawReads.h"
#include "SPGraph.h"
#include "GenomeSegment.h"
#include "ASM.h"
#include "gene.h"
#include "Sim.h"
#include "Acc.h"
#include "SDEAP.h"
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
using namespace std;

int main(int argc, char **argv) {

//run sdeap on BC 

int IDX = atoi(argv[1]);

string gtfdir = "/rhome/ywyang/bigdata/RefSeq/mm9/genes_gtf/";
string tmpdir = "/rhome/ywyang/bigdata/AS/MS/tmp/sdeap/";
string outdir = "/rhome/ywyang/bigdata/AS/MS/pred/sdeap/";

TextIO tio;

vector<string> glistvec = tio.dlread_vector("ms_glist.txt");
string libfname = "ms_sdeap_lib.txt";
vector<string> tmpvec = tio.dlread_vector(libfname);
vector<double> libvec ;
	for(int i =0 ; i< tmpvec.size(); i++){
		libvec.push_back(atof(tmpvec[i].c_str()));
	}
	
string lenfname = "MS_len.txt";	
tmpvec = tio.dlread_vector(lenfname);
vector<int> lenvec ;
for(int i =0 ; i< tmpvec.size(); i++){
	lenvec.push_back(atoi(tmpvec[i].c_str()));
}
		

for(int i =0; i<glistvec.size(); i++){
	if(i%30 == IDX){
		//output bed list
		string gname = glistvec[i];
		vector<string> fnamevec;
		
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E01/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E02/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E03/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E04/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E05/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E06/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E07/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E08/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E09/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E10/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E11/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E12/tophat_out/genes_bed/"+gname);

		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E13/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E14/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E15/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E16/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E17/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E18/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E19/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/E20/tophat_out/genes_bed/"+gname);

		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/M01/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/M02/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/M03/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/M04/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/M05/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/M06/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/M07/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/M08/tophat_out/genes_bed/"+gname);


		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/S01/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/S02/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/S03/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/S04/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/S05/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/S06/tophat_out/genes_bed/"+gname);
		fnamevec.push_back("/rhome/ywyang/bigdata/AS/MS/ES_single_cells/S07/tophat_out/genes_bed/"+gname);


		//run sdeap
		SDEAP sdp;
		sdp.runSDEAP( fnamevec, gname, gtfdir, tmpdir , outdir, libvec, lenvec);
		
	}
}

}