
#include "payoff.hpp"
#include <cmath>
#include <limits>
#include<vector>
#include <algorithm>


//Payoff Class
	//The payoff class is an abstract class
	//The cla ss PayOffCall inherits from payoff, here to price a Call option
	//Regarding parameters, we only need a strike price to define the payoff
	
	//Constructor 
	//Constructor 
PayOffCall::PayOffCall(const double& K)
	:
	 m_K(K)
	{
	};
	
	//Method to compute the initial condition of the Call option
	double PayOffCall::operator() (const double& S,const double& df) const 
	{
		return std::max(S-m_K*df, 0.0); // Call payoff
	};

