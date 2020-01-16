
#include <cmath>
#include <limits>
#include <algorithm>
#include "closed_form.hpp"
//
//Payoff Class
	//PayOff::PayOff(){};
	//We start with the classic call and put payoff
	//Regarding parameters, we only need a strike price to define the payoff
	PayOffCall::PayOffCall(const double& _K):K(_K) {};
	
	double PayOffCall::operator() (const double& S) const 
	{
		return std::max(S-K, 0.0); // Call payoff
	};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mesh Class
	mesh::mesh(const double& spot, const double& maturity,const double& volatility, const long& time_step, const long& steps)
	:dt(maturity/time_step), spot(spot)
	{
		double S0 = log(spot);
		double high_bound = S0 + 5*volatility*sqrt(maturity);
		double low_bound = S0 - 5*volatility*sqrt(maturity);
		
		dx = (high_bound - low_bound)/steps;
		
		//std::vector<double> vector_stock(steps);
		
		for (long i = 0; i < steps+1 ; ++i)
		{
			
			vector_stock.push_back(high_bound + (i - steps)*dx);
		}
		
		//std::vector<double> vector_time(time_step);
		
		for (std::size_t j = 0; j < time_step+1 ; ++j)
		{
			
			vector_time.push_back(j*dt);
			
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
			
	double mesh::getdx() const
	{
		
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

	PDE::PDE(PayOff* _option, const double& dx,const double& dt,const std::vector<double>& time_vector,const std::vector<double>& spot_vector)
	:option(_option)
	{
		size_t nb_step = spot_vector.size();
		//std::vector<double> m_init_vector(nb_step);
		for (std::size_t i = 0; i < nb_step ; ++i)
		{
			m_init_vector.push_back(init_cond(exp(spot_vector[i])));
		}
		
	};

	// Initial condition (vanilla call option), we compute just the payoff created 
	// x parameter stands for the spot
	double PDE::init_cond(const double& x) const 
	{
	  return option->operator()(x);
	};
	std::vector<double> PDE::get_init_vector() const //const forbid to modify the state of my object
    {
		//m_nb_rows = 0;
        return m_init_vector;
    };

    void print(const std::vector<double>& v)
    {
        for(size_t i = 0; i < v.size(); ++i)
        {
            std::cout << v[i] << ",";
        }
        std::cout << std::endl;
    };
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//parameters 

	Parameters::Parameters(double vol, double rate, double theta)
		: pa_vol(vol), pa_Rate(rate), pa_Theta(theta)
	{

	};
	double Parameters::Get_Vol() const{

		return pa_vol;
	};
	double Parameters::Get_Rate() const{

		return pa_Rate;
	};
	double Parameters::Get_Theta() const{

		return pa_Theta;
	};
	Parameters::~Parameters() {

	};
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//bondaries implementation 

	// bound_conditions::bound_conditions()
	// {
		
		// std::cout << "constructor of the bound_conditions" << std::endl;
	// };
	
	//std::vector<double> bound_conditions::operator()(mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman){};
	
	
	
	// std::vector<std::vector<double>>  Neumann::operator()(mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman)
	// {
	
		// double dt = grid.getdt();
		// double dx = grid.getdx(); //need the stock step 
		// double sigma = param.Get_Vol(); //to get the volatility 
		// double rate = param.Get_Rate(); //the get the rate 
		// double theta = param.Get_Theta(); //to get the theta 
		// double size_vec = grid.Getvector_time().size();
		// double maturity = grid.Getvector_time().back();
		// double S0 = grid.get_Spot();
		// double size_spot = grid.Getvector_stock().size();
		
		// double K1 = K_neuman[0];
		// double K2 = K_neuman[1];
		// double K3 = K_neuman[2];
		// double K4 = K_neuman[3];
		
		
		// PDE _payoff(option,dx,dt,grid.Getvector_time(),grid.Getvector_stock()); //create the PDE object from the option 
		
		// std::vector<double>  _init_cond = _payoff.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		
		// double f_0_T = _init_cond[0]; //first element of the vector is the payoff at min S and maturity 
		// double f_N_T = _init_cond.back(); //last element is the payoff at max S and maturity 
		
		// std::vector<double> upper_conditions(size_vec); 
		// std::vector<double> lower_conditions(size_vec); 
		
		// std::fill (upper_conditions.begin(),upper_conditions.end()-1,0);   // we fill the vector with 0 at time 1 to T-1 
		// upper_conditions.back() = f_0_T*exp(-maturity*rate);
		// std::fill (lower_conditions.begin(),lower_conditions.end()-1,0);
		// lower_conditions.back() = f_N_T*exp(-maturity*rate);
		
		// // double beta_left =  dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
		// // double beta_right = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
		
		// // double alpha_left = dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
		// // double alpha = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
		
		// double coef_ =(1 - dt*(1-theta)*rate)/(1+ dt*theta*rate);
		// double coef_K1_K2 = -dt*((-pow(sigma,2)*K1)/2 + (pow(sigma,2)/2 - rate)*K2)/(1+ dt*theta*rate);
		// double coef_K3_K4 = -dt*((-pow(sigma,2)*K3)/2 + ((pow(sigma,2))/2 - rate)*K4)/(1+ dt*theta*rate);
		
		// for (unsigned it = upper_conditions.size(); it != 0; it--)
		// {
		// // reverse iterator to fill the vector from the end to the beginning
			
			// upper_conditions[it]  = coef_*upper_conditions[it-1] + coef_K1_K2; //
			
			// lower_conditions[it] = coef_*lower_conditions[it-1] + coef_K3_K4;
		// }

		// std::vector<std::vector<double>> matrix_neumann(size_spot,std::vector<double> (size_vec));
		  
		// matrix_neumann.front() = upper_conditions;
		
		// for (int i = 1; i < size_spot-1; i++){
			
			// std::vector<double> row_0;
			
			// row_0.resize(size_vec,0.0);
			
			// matrix_neumann.push_back(row_0);
		// }
		
		// matrix_neumann.push_back(lower_conditions);
		
		// std::vector<std::vector<double>> Matrix_conditions(matrix_neumann);
		
		// return matrix_neumann;
	// };
	
	// Derichtlet::Derichtlet(const double& _K): {};
	
	
	//Derichtlet constructor
	// Derichtlet::Derichtlet(const double& _born_sup, const double& _born_min)
	// :
	 // born_sup(_born_sup),
	 // born_min(_born_min)
	// {
	// };
	
	// double PayOffCall::operator() (const double& S) const 
	// {
		// return std::max(S-K, 0.0); // Call payoff
	// };	
	
	bound_conditions::bound_conditions(PDE _payoff, mesh _grid, Parameters _param, PayOff* _option)
	:
	 option(_option),
	 grille(_grid),
	 param(_param),
	 payoff(_payoff)	 
	 {};
	
	Derichtlet::Derichtlet(PDE payoff, mesh grid, Parameters param, PayOff* option)
	:bound_conditions(payoff,grille,param,option)
	{
		
		double dt = grid.getdt();
		double dx = grid.getdx(); //need the stock step 
		double sigma = param.Get_Vol(); //to get the volatility 
		double rate = param.Get_Rate(); //the get the rate 
		double theta = param.Get_Theta(); //to get the theta 
		size_t size_vec = grid.Getvector_time().size();
		double maturity = grid.Getvector_time().back();
		double S0 = grid.get_Spot();
		size_t size_spot = grid.Getvector_stock().size();
		double r = param.Get_Rate();
		std::vector<double>  _init_cond = payoff.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		
		//double f_0_T = _init_cond[0]*exp(-maturity*rate); //first element of the vector is the payoff at min S and maturity 
		double f_0_T = _init_cond[0];
		// double f_N_T = _init_cond.back()*exp(-maturity*rate); //last element is the payoff at max S and maturity 
		double f_N_T = _init_cond.back();
		
		matrix_derichtlet.resize(size_vec*2); 
		//std::vector<double> lower_conditions(size_vec); 
		
		for(size_t i = 0; i < size_vec; ++i)
        {
            matrix_derichtlet[i] = f_0_T*exp(-r*dt*i); // for here
			matrix_derichtlet[size_vec+i] = f_N_T*exp(-r*dt*i);
        }
		 //matrix_derichtlet = upper_conditions;
		// matrix_derichtlet.front() = upper_conditions;
		
		// for (int i = 1; i < size_spot-1; i++)
		// {
			
			// std::vector<double> row_0;
			
			// row_0.resize(size_vec,0.0);
			
			// matrix_derichtlet.push_back(row_0);
		// }
		
		// matrix_derichtlet.push_back(lower_conditions);
		
		//std::vector<std::vector<double>> Matrix_conditions(matrix_derichtlet);

		//return matrix_derichtlet;
		
			
	};
	
	std::vector<double> Derichtlet::get_cond() const
	{
		return matrix_derichtlet;
	};
	
	Neumann::Neumann(PDE payoff, mesh grid, Parameters param, PayOff* option,std::vector<double>& const_vector)
	:bound_conditions(payoff,grille,param,option)
	{
		
		double dt = grid.getdt();
		double dx = grid.getdx(); //need the stock step 
		double sigma = param.Get_Vol(); //to get the volatility 
		double rate = param.Get_Rate(); //the get the rate 
		double theta = param.Get_Theta(); //to get the theta 
		size_t size_vec = grid.Getvector_time().size();
		double maturity = grid.Getvector_time().back();
		double S0 = grid.get_Spot();
		size_t size_spot = grid.Getvector_stock().size();
		double r = param.Get_Rate();
		std::vector<double>  _init_cond = payoff.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		
		//double f_0_T = _init_cond[0]*exp(-maturity*rate); //first element of the vector is the payoff at min S and maturity 
		double f_0_T = _init_cond[0];
		// double f_N_T = _init_cond.back()*exp(-maturity*rate); //last element is the payoff at max S and maturity 
		double f_N_T = _init_cond.back();
		
		   double K1 = const_vector[0];
		   double K2 = const_vector[1];
		   double K3 = const_vector[2];
		   double K4 = const_vector[3];
		   
		double coef_ =(1 - dt*(1-theta)*rate)/(1+ dt*theta*rate);
		double coef_K1_K2 = -dt*((-pow(sigma,2)*K1)/2 + (pow(sigma,2)/2 - rate)*K2)/(1+ dt*theta*rate);
		double coef_K3_K4 = -dt*((-pow(sigma,2)*K3)/2 + (pow(sigma,2)/2 - rate)*K4)/(1+ dt*theta*rate);
		
		matrix_neumann.resize(size_vec);
		matrix_neumann.push_back(f_0_T*exp(-maturity*rate));
		matrix_neumann.resize(size_vec*2);
		matrix_neumann.push_back(f_N_T*exp(-maturity*rate));
		
		for (size_t it = size_vec-1; it != 0; it--)
		{
		      //reverse iterator to fill the vector from the end to the beginning
			  
			  
			  size_t tt = size_vec + it;
			  matrix_neumann[it]  = (coef_*matrix_neumann[it-1] + coef_K1_K2)*exp(-r*dt*it);
			
			  matrix_neumann[tt] = (coef_*matrix_neumann[tt-1] + coef_K3_K4)*exp(-r*dt*it);
		 }

		// std::vector<double> upper_conditions(size_vec); 
		// std::vector<double> lower_conditions(size_vec); 
		
		// std::fill (upper_conditions.begin(),upper_conditions.end()-1,0);   // we fill the vector with 0 at time 1 to T-1 
		// upper_conditions.back() = f_0_T*exp(-maturity*rate);
		// std::fill (lower_conditions.begin(),lower_conditions.end()-1,0);
		// lower_conditions.back() = f_N_T*exp(-maturity*rate);
		
		// // double beta_left =  dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
		// // double beta_right = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
		
		// // double alpha_left = dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
		// // double alpha = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
		

		
		// for (unsigned it = upper_conditions.size(); it != 0; it--)
		// {
		// // reverse iterator to fill the vector from the end to the beginning
			
			// upper_conditions[it]  = coef_*upper_conditions[it-1] + coef_K1_K2; //
			
			// lower_conditions[it] = coef_*lower_conditions[it-1] + coef_K3_K4;
		// }

		// std::vector<std::vector<double>> matrix_neumann(size_spot,std::vector<double> (size_vec));
		  
		// matrix_neumann.front() = upper_conditions;
		
		// for (int i = 1; i < size_spot-1; i++){
			
			// std::vector<double> row_0;
			
			// row_0.resize(size_vec,0.0);
			
			// matrix_neumann.push_back(row_0);
		// }
		
		// matrix_neumann.push_back(lower_conditions);
		
		// std::vector<std::vector<double>> Matrix_conditions(matrix_neumann);
		
		// return matrix_neumann;
		
			
	};
	
    std::vector<double> Neumann::get_cond() const
	{
		return matrix_neumann;
	};
	
	// std::vector<std::vector<double>>   bound_conditions::boundaries_compute(mesh grid, Parameters param, PayOff* option, bound_conditions* bound_func, std::vector<double> K_neuman){
		
	// return (*bound_func)(grid, param, option, K_neuman);
		
	// };
	
//Autre version code Dirichlet///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	
	// std::vector<std::vector<double>>  Derichtlet::operator()(PDE _payoff, mesh grid, Parameters param, PayOff* option)
	// {
		
		// double dt = grid.getdt();
		// double dx = grid.getdx(); //need the stock step 
		// double sigma = param.Get_Vol(); //to get the volatility 
		// double rate = param.Get_Rate(); //the get the rate 
		// double theta = param.Get_Theta(); //to get the theta 
		// size_t size_vec = grid.Getvector_time().size();
		// double maturity = grid.Getvector_time().back();
		// double S0 = grid.get_Spot();
		// size_t size_spot = grid.Getvector_stock().size();
		
		// //PDE _payoff(option,dx,dt,grid.Getvector_time(),grid.Getvector_stock()); //create the PDE object from the option 
		
		// std::vector<double>  _init_cond = _payoff.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		
		// //double f_0_T = _init_cond[0]*exp(-maturity*rate); //first element of the vector is the payoff at min S and maturity 
		// double f_0_T = _init_cond[0];
		// // double f_N_T = _init_cond.back()*exp(-maturity*rate); //last element is the payoff at max S and maturity 
		// double f_N_T = _init_cond.back();
		
		// std::vector<double> upper_conditions(size_vec); 
		// std::vector<double> lower_conditions(size_vec); 
		
		// std::fill (upper_conditions.begin(),upper_conditions.end(),f_0_T);   // we fill the vector with the terminal condition as Vn f_n is a constant (l or h)  
		// std::fill (lower_conditions.begin(),lower_conditions.end(),f_N_T);
		
		// //double beta_left =  dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
		// //double beta_right = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
		
		// //double alpha_left = dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
		// //double alpha = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
		
		// // double coef_ =1;
		// // double coef_K1_K2 = 0;
		// // double coef_K3_K4 = 0;
		
		// // for (unsigned it = upper_conditions.size(); it != 0; it--){
		// // //reverse iterator to fill the vector from the end to the beginning
			
			// // upper_conditions[it]  = coef_*upper_conditions[it-1] + coef_K1_K2; //
			
			// // lower_conditions[it] = coef_*lower_conditions[it-1] + coef_K3_K4;
		// // }
		
		// std::vector<std::vector<double>> matrix_derichtlet(size_spot,std::vector<double> (size_vec));
		  
		// matrix_derichtlet.front() = upper_conditions;
		
		// for (int i = 1; i < size_spot-1; i++){
			
			// std::vector<double> row_0;
			
			// row_0.resize(size_vec,0.0);
			
			// matrix_derichtlet.push_back(row_0);
		// };
		
		// matrix_derichtlet.push_back(lower_conditions);
		
		// //std::vector<std::vector<double>> Matrix_conditions(matrix_derichtlet);

		// return matrix_derichtlet;
		
			
	// };
	
	// std::vector<std::vector<double>>   bound_conditions::boundaries_compute(mesh grid, Parameters param, PayOff* option, bound_conditions* bound_func, std::vector<double> K_neuman){
		
	// return (*bound_func)(grid, param, option, K_neuman);
		
	// };
	
	
	// const std::vector<std::vector<double>> bound_conditions::get_matrix(){
		
		// return Matrix_conditions;
	// };
	
	
	// bound_conditions::~bound_conditions(){}; //destructor of the boundaries_conditions 
	


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
		