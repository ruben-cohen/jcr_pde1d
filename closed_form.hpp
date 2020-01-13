#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP
#include <iostream>
#include <vector>
#include <algorithm> // Act on containers through iterators to apply modyfing/non_modifying operations



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
	
	class PDE
	{
		public:
		PDE(PayOff* _option);
		PayOff* option;
		
		// std::vector<double> Mid_diag_coeff(mesh_spot grid, parameters param,bool A); 
		// std::vector<double> Upper_diag_coeff(mesh_spot grid, parameters param,bool A);
		// std::vector<double> Lower_diag_coeff(mesh_spot grid, parameters param,bool A);

		//all in function to get the price of the option and the greeks ? 
		//std::vector<double> resolution();
		
		//function to compute the greeks ? 

		//double boundary_left(double t, double x) const;
		//double boundary_right(double t, double x) const;
		
		double init_cond(const double& x) const; //To compute the payoff we want to implement (at maturity)
		//std::vector<double> init_cond(
		
		//std::vector<double> CranckNicholson_algo(mesh_spot grid, parameters param); //main function where we define the CN procedure
		//std::vector<double> Thomas_triLinMatrix_inverse(); //need to inverse A in the system that is a diagonal matrix 
		
		//set of function to define the A matrix (3 better than just one huge matrix ?) 
		//in each function we only need to input the grid and the parameters as we will get theta, sigma and r from parameters class 
		// we will get dx and dt from grid class 
	};
	
// class VanillaOption 
	// {
		 // public:
		  // PayOff* pay_off; //Pointer 

		  // double K;
		  // double r;
		  // double T;
		  // double sigma;

		  // VanillaOption();
		  // VanillaOption(double _K, double _r, double _T, 
						// double _sigma, PayOff* _pay_off);
	// };
	
	    // double vanilla_payoff(double fwd, double strike, bool is_call);
    // double bs_time_value(double fwd, double strike, double volatility, double maturity);
    // double bs_price(double fwd, double strike, double volatility, double maturity, bool is_call);

    // std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call);
    // std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity);
    // std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, bool is_call);

#endif
