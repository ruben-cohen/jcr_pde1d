#ifndef PAY_OFF_CPP
#define PAY_OFF_CPP

#include "payoff.hpp"

namespace projet
{
	
	PayOff::PayOff() 
	{	
	}
	
//We start with the classic call and putpayoff
	//Regarding parameters, we only need a strike price to define de payoff of both options
	PayOffCall::PayOffCall(const double& _K) 
	{
		K = _K; 
	}
	double PayOffCall::operator() (const double& S) const 
	{
		return std::max(S-K, 0.0); // Call payoff
	}
	PayOffPut::PayOffPut(const double& _K)
	{
		K = _K;
	}

	// Over-ridden operator() method, which turns PayOffPut into a function object
	double PayOffPut::operator() (const double& S) const
	{
	  return std::max(K-S, 0.0); // Put pay-off
	}

}
#endif