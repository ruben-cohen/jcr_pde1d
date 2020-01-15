#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP
#include <iostream>
#include <vector>
#include <algorithm> // Act on containers through iterators to apply modyfing/non_modifying operations
//

//Payoff Class
	class PayOff 
	{
		 public:
		  PayOff(){};
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
		  PayOffCall(const double& _K);
		  virtual ~PayOffCall() {};

		  // Virtual function is now over-ridden (not pure-virtual anymore)
		  virtual double operator() (const double& S) const;
		  //virtual double init_cond(const double& S)) const;
		  
		 private:
		  double K; // Variable for the Strike price

		 
	};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mesh Class	
	class mesh
	{
	
	public: 
	mesh(const double& spot, const double& maturity, const double& volatility,const long& time_step,const size_t& steps);
	~mesh();
	std::vector<double> Getvector_time() const;
	std::vector<double> Getvector_stock() const;
	double getdx() const;
	double getdt() const;
	double get_Spot() const;
	
	
	private: 
	
	std::vector<double> vector_time;
	std::vector<double> vector_stock;
	double dx;
	double dt;
	double spot;
	};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// PDE Solver	
	class PDE
	{
		public: 
		
		PDE(PayOff* _option, const double& dx,const double& dt,const std::vector<double>& time_vector,const std::vector<double>& spot_vector);
		PayOff* option;
		
		//To compute the payoff we want to implement (at maturity)
		double init_cond(const double& x) const;
		std::vector<double> get_init_vector() const;
		void print(const std::vector<double>& v);
		// std::vector<double> init_cond(		
		// std::vector<double> Mid_diag_coeff(mesh_spot grid, parameters param,bool A); 
		// std::vector<double> Upper_diag_coeff(mesh_spot grid, parameters param,bool A);
		// std::vector<double> Lower_diag_coeff(mesh_spot grid, parameters param,bool A);

		// all in function to get the price of the option and the greeks ? 
		// std::vector<double> resolution();
		
		// function to compute the greeks ? 

		// double boundary_left(double t, double x) const;
		// double boundary_right(double t, double x) const;
		
		private:
		
		std::vector<double> m_init_vector;
		
		
		
		// std::vector<double> CranckNicholson_algo(mesh_spot grid, parameters param); //main function where we define the CN procedure
		// std::vector<double> Thomas_triLinMatrix_inverse(); //need to inverse A in the system that is a diagonal matrix 
		
		// set of function to define the A matrix (3 better than just one huge matrix ?) 
		// in each function we only need to input the grid and the parameters as we will get theta, sigma and r from parameters class 
		// we will get dx and dt from grid class 
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
	
	class bound_conditions {
	public: // all virtuals ? and we have all through neumann et/ou derichtlet ? 
	
	bound_conditions();
	
	std::vector<std::vector<double>>  operator() (mesh grid, Parameters param, PayOff* option, std::vector<double> K_neuman);
	~bound_conditions();
	
	const std::vector<std::vector<double>> get_matrix();
	
	static std::vector<std::vector<double>>   boundaries_compute(mesh grid, Parameters param, PayOff* option, bound_conditions* bound_func, std::vector<double> K_neuman ={0,0,0,0});
	//this function takes the same parameters as the bound_conditions + the type of boundaries we want to compute and will return the appropriate boundaries.

	private:
	
	std::vector<std::vector<double>> Matrix_conditions;
		
	};
	
	class Neumann : public  bound_conditions {
		
	public:
	Neumann(){};
		
	std::vector<std::vector<double>>  operator() (mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);
	
	
	private: 
	
	std::vector<double> K_neuman; 
	
	};
	
	class Derichtlet: public bound_conditions {
		
	public:
	
	Derichtlet(){};
	
	std::vector<std::vector<double>>  operator() (mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);

	private: 
	
	std::vector<double> K_neuman; 
		
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
