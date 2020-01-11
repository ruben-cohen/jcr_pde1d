#ifndef VANILLA_OPTION_H
#define VANILLA_OPTION_H

#include "payoff.hpp"


	class VanillaOption 
	{
		 public:
		  PayOff* pay_off; //Pointer 

		  double K;
		  double r;
		  double T;
		  double sigma;

		  VanillaOption();
		  VanillaOption(double _K, double _r, double _T, 
						double _sigma, PayOff* _pay_off);
	};


#endif