
#include <cmath>
#include <limits>
#include <algorithm>

#include "closed_form.hpp"

//Payoff Class
	PayOff::PayOff() 
	{	
	}
	//We start with the classic call and put payoff
	//Regarding parameters, we only need a strike price to define de payoff of both options
	PayOffCall::PayOffCall(const double& _K) 
	{
		K = _K; 
	}
	double PayOffCall::operator() (const double& S) const 
	{
		return std::max(S-K, 0.0); // Call payoff
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mesh Class
	mesh::mesh(const double& spot, const double& maturity,const double& volatility, const long& time_step, const size_t& steps)
	:dt(time_step), spot(spot) //dx(steps)
	{
		double S0 = log(spot);
		double high_bound = S0 + 5*volatility*sqrt(maturity);
		double low_bound = S0 - 5*volatility*sqrt(maturity);
		
		double dx_spot = (high_bound - low_bound)/steps;
		
		std::vector<double> vector_stock(steps);
		
		for (std::size_t i = 0; i < steps ; ++i){
			
			vector_stock[i] = S0 + (i - steps)* dx_spot;
		}
		
		std::vector<double> vector_time(time_step);
		
		for (std::size_t j = 0; j < maturity ; ++j){
			
			vector_time[j] = j*time_step;
			
		dx = dx_spot;
		}
	};
	mesh::~mesh() {};
			//std::cout << "destructor of the grid" << std::endl; //destructor of the grid 
	
	std::vector<double> mesh::Getvector_time()const{
		
			return vector_time;
	}; //useful to get the vector of time from the mesh 
			
	std::vector<double> mesh::Getvector_stock() const{
		
		return vector_stock;
	
	}; //useful to get the vector of stock path from the mesh 
			
	double mesh::getdx() const{
		
		return dx;
	
	}; //we will need to get the dx for CN algo 
	
	double mesh::getdt() const{
		
		return dt;
	
	}; // we will need the dt for CN algo 
	
	double mesh::get_Spot() const{
		
		return spot;
	};
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// PDE Solver

	PDE::PDE(PayOff* _option)
	:option(_option)
	{};

	// Initial condition (vanilla call option), we compute just the payoff created 
	// x parameter stands for the spot
	double PDE::init_cond(const double& x) const 
	{
	  return option->operator()(x);
	}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Poubelle
	

    // double ncdf(double x)
    // {
        // return 0.5 * std::erfc(-x / std::sqrt(2));
    // }
    
    // double vanilla_payoff(double fwd, double strike, bool is_call)
    // {
        // return std::max(is_call ? fwd - strike : strike - fwd, 0.);
    // }

    // double bs_time_value(double fwd, double strike, double volatility, double maturity)
    // {
        // if(strike == 0.)
        // {
            // return 0.;
        // }
        // else
        // {
            // double stddev = volatility * std::sqrt(maturity);
            // if(stddev == 0.)
            // {
                // return 0.;
            // }
            // double tmp = std::log(fwd / strike) / stddev;
            // double d1 = tmp + 0.5 * stddev;
            // double d2 = tmp - 0.5 * stddev;
            // double res;
            // if(fwd > strike)
            // {
                // res = strike * ncdf(-d2) - fwd * ncdf(-d1);
            // }
            // else
            // {
                // res = fwd * ncdf(d1) - strike * ncdf(d2);
            // }
            // if(res <= std::numeric_limits<double>::min())
            // {
                // res = 0.;
            // }
            // return res;
        // }
    // }

    // double bs_price(double fwd, double strike, double volatility, double maturity, bool is_call)
    // {
        // return vanilla_payoff(fwd, strike, is_call) + bs_time_value(fwd, strike, volatility, maturity);
    // }

    // std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call)
    // {
        // std::vector<double> res(fwd.size());
        // for(std::size_t i = 0; i < fwd.size(); ++i)
        // {
            // res[i] = vanilla_payoff(fwd[i], strike, is_call);
        // }
        // return res;
    // }

    // std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity)
    // {
        // std::vector<double> res(fwd.size());
        // for(std::size_t i = 0; i < fwd.size(); ++i)
        // {
            // res[i] = bs_time_value(fwd[i], strike, volatility, maturity);
        // }
        // return res;
    // }

    // std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, bool is_call)
    // {
        // std::vector<double> res(fwd.size());
        // for(std::size_t i = 0; i < fwd.size(); ++i)
        // {
            // res[i] = bs_price(fwd[i], strike, volatility, maturity, is_call);
        // }
        // return res;
    // }
	// VanillaOption::VanillaOption()
	 // {
	 // }

	// VanillaOption::VanillaOption(double _K, double _r, double _T, double _sigma, PayOff* _pay_off)
		// : 
		// K(_K), r(_r), T(_T), sigma(_sigma), pay_off(_pay_off)
		// {
		// }
		