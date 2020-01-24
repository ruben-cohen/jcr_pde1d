#ifndef PAY_OFF_CPP
#define PAY_OFF_CPP

#include "payoff.hpp"

namespace projet
{
	
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
}
#endif
