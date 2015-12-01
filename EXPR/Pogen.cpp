
#include "Pogen.h"
using namespace boost;

double  Pogen::NB_Gen(int n, double p)
{
	static boost::mt19937 randgen(0);
	boost::mt19937 eng;
	boost::random::negative_binomial_distribution<> pd(n ,p);
	static boost::variate_generator <boost::mt19937&,boost::random::negative_binomial_distribution <> > generator(randgen, pd);

  // sample from the distribution
  return (double) generator();
  
  
}

double  Pogen::gen(double mean , double dispersion){
	 
	//double variance = mean*mean*dispersion;
	
	//double p = mean/variance;
  double size = 1/dispersion;
  double prob =  prob = size/(size+mean);
  //std::cout << "n,p = " << n<< ","<< p<< ","<< mean<< std::endl;
	return NB_Gen(size,prob);
	
}


