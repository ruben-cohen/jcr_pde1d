#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP
#include <iostream>
#include <vector>
#include <algorithm> // Act on containers through iterators to apply modyfing/non_modifying operations
//Allows to create and manipulate xarrays
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xeval.hpp"
//Allows to perform linear algebra operations on xarrays
#include "xtensor-blas/xlinalg.hpp"
////

//Payoff Class
	class PayOff 
	{
		 public:
			//Constructor
			PayOff(){};
			// Virtual destructor to avoid memory leaks when destroying the base and inherited classes
			virtual ~PayOff() {}; 
			// We turn the payoff into a functor, goal is to compute initial condition
			virtual double operator() (const double& S) const = 0;
			//virtual double init_cond(const double& S)) const;
	};
	
	
	
	//First payoff class created for european call options
	//It inherits from Payoff
	class PayOffCall : public PayOff 
	{
		public:
		
			PayOffCall(const double& K);
			virtual ~PayOffCall() {};

			// Virtual function is now over-ridden (not pure-virtual anymore)
			virtual double operator() (const double& S) const;
			//virtual double init_cond(const double& S)) const;
		  
		 private:
		 
			double m_K; // Variable for the Strike price

		 
	};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mesh Class	
	class mesh
	{
	
	public: 
	
		mesh(const double& spot, const double& maturity, const double& volatility,const long& time_step,const long& steps, PayOff* _option);
		~mesh();
		PayOff* option;
		
		std::vector<double> Getvector_time() const;
		std::vector<double> Getvector_stock() const;
		double getdx() const;
		double getdt() const;
		double get_Spot() const;
		double init_cond(const double& x) const;
		std::vector<double> get_init_vector() const;
	
	
	private: 
	
		std::vector<double> m_vector_time;
		std::vector<double> m_vector_stock;
		double m_dx;
		double m_dt;
		double m_spot;
		long m_steps;
		long m_time_step;
		std::vector<double> m_init_vector;
		
	};
	
//Method to print vector content
void print(const std::vector<double>& v);
//Method to transform the vector boundaries into matrix for resolution;
std::vector<double> transform_matrix(std::vector<double> vector_init, double nb_rows);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// paramaters 
	
	class Parameters { 
	public:
		Parameters(double vol, double rate, double theta);
		double Get_Vol() const;
		double Get_Rate() const;
		double Get_Theta() const;
		~Parameters();

	private:
		
		double pa_vol;
		double pa_Rate;
		double pa_Theta;
	};
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// boundaries class
	
	class bound_conditions 
	{
		public:
		
			bound_conditions(mesh _grid, Parameters _param);
			
			//std::vector<double> cond() const;
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 //std::vector<double>  cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
		
		protected:
		
			mesh m_grille;
			Parameters m_param;
	
		//private:
		
		//std::vector<std::vector<double>> Matrix_conditions;
			
	};
	class Derichtlet: public bound_conditions 
	{
		
		public:
		
			Derichtlet(mesh _grid, Parameters _param);
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);
			 //std::vector<double> cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 std::vector<double> get_cond() const;

			std::vector<double> matrix_derichtlet;
			
		//std::vector<double> K_neuman; 
		
	};
	
	
	class Neumann: public bound_conditions 
	{
		
		public:
		
			Neumann(mesh _grid, Parameters _param, std::vector<double>& const_vector);
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);
			 //std::vector<double> cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 std::vector<double> get_cond() const;

			std::vector<double> matrix_neumann;
			
		//std::vector<double> K_neuman; 
		
	};
	
class Solve
    {
        
    public:

      Solve(mesh _grid, Parameters _param, std::vector<double>& conditions);
		
		~Solve();
		
		xt::xarray<double> Get_FX_n() const; //to get the price at step n;
		
	protected:
		
	mesh m_grille;
	Parameters m_param;
		
	private: 
	
	xt::xarray<double> _FX_n; // vector of price at n 
	std::vector<double> conditions;
        
        //double dt; //to get the time step
        //double dx; //to get the stock value step
        //int nb_step_time; //to get the number of time steps on which we iterate
        //double v; //to get the volatility
        //double r; //to get the rate
        //double theta; //to get the theta
        //int nb_step_spot; //to get the size of the matrix A and B (N-1)


        //Iterates on each time step of the grid to solve the equation. Beginning from time T-1 and solving for time 1

    };
	
	
	



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Poubelle

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// // PDE Solver	
	// // class PDE
	// // {
		// // public: 
			
			// // PDE(PayOff* _option, const double& dx,const double& dt,const std::vector<double>& time_vector,const std::vector<double>& spot_vector);
			// // PayOff* option;
			// // ~PDE();
			// //To compute the payoff we want to implement (at maturity)
			// // double init_cond(const double& x) const;
			// // std::vector<double> get_init_vector() const;
			
		// // private:
		
			// // std::vector<double> m_init_vector;
		
	
	// // };

	
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
