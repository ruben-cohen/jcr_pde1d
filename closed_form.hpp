#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP
#include <iostream>
#include <vector>
#include <algorithm> // Act on containers through iterators to apply modyfing/non_modifying operations

    double vanilla_payoff(double fwd, double strike, bool is_call);
    double bs_time_value(double fwd, double strike, double volatility, double maturity);
    double bs_price(double fwd, double strike, double volatility, double maturity, bool is_call);

    std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call);
    std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity);
    std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, bool is_call);

	class PayOff 
	{
		 public:
		  PayOff();
		  // Virtual destructor to avoid memory leaks when destroying the base and inherited classes
		  virtual ~PayOff() {}; 
		  // We turn the class into a functor (object we can call just like an object)
		  virtual double operator() (const double& S) const = 0;
		  //virtual double init_cond(const double& S)) const;
	};
	//All the following class are bound to the Payoff class
	// first payoff class created for european call options
	class PayOffCall : public PayOff 
	{
		public:
		  PayOffCall(const double& K_);
		  virtual ~PayOffCall() {};

		  // Virtual function is now over-ridden (not pure-virtual anymore)
		  virtual double operator() (const double& S) const;
		  //virtual double init_cond(const double& S)) const;
		  
		 private:
		  double K; // Variable for the Strike price

		 
	};
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
