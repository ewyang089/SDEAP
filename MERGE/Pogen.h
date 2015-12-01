#ifndef POGEN_H
#define POGEN_H

#include <math.h>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random.hpp>
#include <boost/random/negative_binomial_distribution.hpp>

class Pogen{
	public:
	
	double  NB_Gen(int n, double p);
	double  gen(double mean ,double variance );
			
		
};

#endif
