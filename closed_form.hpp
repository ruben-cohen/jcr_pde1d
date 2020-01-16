#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP
#include <iostream>
#include <vector>
#include <algorithm> // Act on containers through iterators to apply modyfing/non_modifying operations
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
	
		mesh(const double& spot, const double& maturity, const double& volatility,const long& time_step,const long& steps);
		~mesh();
		
		std::vector<double> Getvector_time() const;
		std::vector<double> Getvector_stock() const;
		double getdx() const;
		double getdt() const;
		double get_Spot() const;
	
	
	private: 
	
		std::vector<double> m_vector_time;
		std::vector<double> m_vector_stock;
		double m_dx;
		double m_dt;
		double m_spot;
		long m_steps;
		long m_time_step;
	};
	
void print(const std::vector<double>& v);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// PDE Solver	
	class PDE
	{
		public: 
			
			PDE(PayOff* _option, const double& dx,const double& dt,const std::vector<double>& time_vector,const std::vector<double>& spot_vector);
			PayOff* option;
			~PDE();
			//To compute the payoff we want to implement (at maturity)
			double init_cond(const double& x) const;
			std::vector<double> get_init_vector() const;
			
		private:
		
			std::vector<double> m_init_vector;
		
	
	};

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
		
			bound_conditions(PDE _payoff, mesh _grid, Parameters _param, PayOff* _option);
			
			//std::vector<double> cond() const;
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 //std::vector<double>  cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
		
		protected:
		
			PDE m_payoff;
			mesh m_grille;
			Parameters m_param;
			PayOff* m_option;



		//private:
		
		//std::vector<std::vector<double>> Matrix_conditions;
			
	};
	class Derichtlet: public bound_conditions 
	{
		
		public:
		
			Derichtlet(PDE _payoff, mesh _grid, Parameters _param, PayOff* _option);
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);
			 //std::vector<double> cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 std::vector<double> get_cond() const;

			std::vector<double> matrix_derichtlet;
			
		//std::vector<double> K_neuman; 
		
	};
	
	
	class Neumann: public bound_conditions 
	{
		
		public:
		
			Neumann(PDE _payoff, mesh _grid, Parameters _param, PayOff* _option, std::vector<double>& const_vector);
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);
			 //std::vector<double> cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 std::vector<double> get_cond() const;

			std::vector<double> matrix_neumann;
			
		//std::vector<double> K_neuman; 
		
	};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Poubelle
	
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
