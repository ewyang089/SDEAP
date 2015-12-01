#ifndef CHRM_H
#define CHRM_H

#include "GenomeSegment.h"


class Chrm{
	public:
	GenomeSegment gs; 
	void initChrm(std::string fname, int chr_len);
	void mapReads_BED(std::string fname);
	void printReadCounts(std::string outfname);
	private:
	int	_check_split(std::string CICAR);
		
		
};

#endif