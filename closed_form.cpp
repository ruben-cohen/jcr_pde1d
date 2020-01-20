
#include <cmath>
#include <limits>
#include <algorithm>
#include "closed_form.hpp"


//Payoff Class
	//The payoff class is an abstract class
	//The class PayOffCall inherits from payoff, here to price a Call option
	//Regarding parameters, we only need a strike price to define the payoff
	
	//Constructor 
	PayOffCall::PayOffCall(const double& K)
	:
	 m_K(K)
	{
	};
	
	//Method to compute the initial condition of the Call option
	double PayOffCall::operator() (const double& S) const 
	{
		return std::max(S-m_K, 0.0); // Call payoff
	};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mesh Class
	//This class builds the mesh, giving dt, dS, and the axes of the mesh as output
	mesh::mesh(const double& spot, const double& maturity,const double& volatility, const long& time_step, const long& steps,PayOff* _option)
	:
	 option(_option),
	 m_dt(maturity/time_step),
	 m_spot(spot),
	 m_steps(steps),
	 m_time_step(time_step)
	{
		//The mesh is centered around log(S0)
		double S0 = log(m_spot);
		double high_bound = S0 + 5*volatility*sqrt(maturity);
		double low_bound = S0 - 5*volatility*sqrt(maturity);
		
		//dS
		m_dx = (high_bound - low_bound)/steps;
		
		//Axe of spot price
		for (long i = 0; i < m_steps+1 ; ++i)
		{
			
			m_vector_stock.push_back(high_bound + (i - m_steps)*m_dx);
		}
		//Axe of time
		for (std::size_t j = 0; j < m_time_step+1 ; ++j)
		{
			m_vector_time.push_back(j*m_dt);	
		}
		
		size_t nb_step = m_vector_stock.size();
		//std::vector<double> m_init_vector(nb_step);
		
		for (std::size_t i = 0; i < nb_step ; ++i)
		{
			m_init_vector.push_back(init_cond(exp(m_vector_stock[i])));
		}
	};

		
	//Multiple Get methods to return private variables of the class


	//useful to get the vector of time from the mesh 
	std::vector<double> mesh::Getvector_time()const
	{
		return m_vector_time;
	}; 
	//useful to get the vector of stock path from the mesh 
	std::vector<double> mesh::Getvector_stock() const
	{
		return m_vector_stock;
	}; 
	double mesh::init_cond(const double& x) const 
	{
	  return option->operator()(x);
	};
	
	//we will need to get the dx for CN algo	
	double mesh::getdx() const
	{
		return m_dx;
	}; 
	std::vector<double> mesh::get_init_vector() const //const forbid to modify the state of my object
    {
        return m_init_vector;
    };
	
	// we will need the dt for CN algo 
	double mesh::getdt() const
	{
		return m_dt;
	}; 
	
	double mesh::get_Spot() const
	{
		return m_spot;
	};
	
	//Destructor
	mesh::~mesh() {};
	
    void print(const std::vector<double>& v)
    {
        for(size_t i = 0; i < v.size(); ++i)
        {
            std::cout << v[i] << ","; 
        }
        std::cout << std::endl;
    };

	std::vector<double> transform_matrix(std::vector<double> vector_init, double nb_rows){
	
	double endv; 
	std::vector<double> upper_bound(vector_init);
	std::vector<double> lower_bound(vector_init);
	 
	 if (vector_init.size()%2==0){endv=vector_init.size()%2;}
	 else {endv=std::floor(vector_init.size()/2) -1;};
	 
	 upper_bound.resize(endv);
	 lower_bound.resize(-endv -1);
	 
	std::vector<double> matrix(upper_bound);
	
	for(int i=0; i<endv*nb_rows; i++){
		
		matrix.push_back(0.0);
	};
	
	for(size_t t=0; t<lower_bound.size(); t++)
	{
		matrix.push_back(lower_bound[t]);
	};
	 
	return matrix;		
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//parameters 

	Parameters::Parameters(double vol, double rate, double theta)
		:
		 pa_vol(vol), 
		 pa_Rate(rate), 
		 pa_Theta(theta)
	{
	};
	double Parameters::Get_Vol() const
	{
		return pa_vol;
	};
	double Parameters::Get_Rate() const
	{
		return pa_Rate;
	};
	double Parameters::Get_Theta() const
	{
		return pa_Theta;
	};
	Parameters::~Parameters() {

	};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Boundaries	
	
	bound_conditions::bound_conditions(mesh _grid, Parameters _param)
	:
	 m_grille(_grid),
	 m_param(_param)
	 {};
	 
	
	Derichtlet::Derichtlet(mesh m_grid, Parameters m_param)
	:
	 bound_conditions(m_grid,m_param)
	{
		//From mesh object
		double dt = m_grid.getdt();
		double dx = m_grid.getdx(); //need the stock step 
		size_t size_vec = m_grid.Getvector_time().size();
		double maturity = m_grid.Getvector_time().back();
		double S0 = m_grid.get_Spot();
		size_t size_spot = m_grid.Getvector_stock().size();
		
		//From parametre object
		double sigma = m_param.Get_Vol(); //to get the volatility 
		double rate = m_param.Get_Rate(); //the get the rate 
		double theta = m_param.Get_Theta(); //to get the theta 
		double r = m_param.Get_Rate();
		
		//From PDE object
		std::vector<double>  _init_cond = m_grid.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		double f_0_T = _init_cond[0];
		double f_N_T = _init_cond.back();
		
		matrix_derichtlet.resize(size_vec*2); 
		//std::vector<double> lower_conditions(size_vec); 
		
		for(size_t i = 0; i < size_vec; ++i)
        {
            matrix_derichtlet[i] = f_0_T*exp(-r*dt*i); // for here
			matrix_derichtlet[size_vec+i] = f_N_T*exp(-r*dt*i);
        }
	};
	
	std::vector<double> Derichtlet::get_cond() const
	{
		return matrix_derichtlet;
	};
	
	Neumann::Neumann( mesh m_grid, Parameters m_param,std::vector<double>& const_vector)
	:bound_conditions(m_grille,m_param)
	{
		//From mesh object
		double dt = m_grid.getdt();
		double dx = m_grid.getdx(); //need the stock step 
		size_t size_vec = m_grid.Getvector_time().size();
		double maturity = m_grid.Getvector_time().back();
		double S0 = m_grid.get_Spot();
		size_t size_spot = m_grid.Getvector_stock().size();		
		
		//From parametre object
		double sigma = m_param.Get_Vol(); //to get the volatility 
		double rate = m_param.Get_Rate(); //the get the rate 
		double theta = m_param.Get_Theta(); //to get the theta 
		double r = m_param.Get_Rate();
		
		//From PDE object
		std::vector<double>  _init_cond = m_grid.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		double f_0_T = _init_cond[0];
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
	};
	
    std::vector<double> Neumann::get_cond() const
	{
		return matrix_neumann;
	};
	
	
Solve::Solve(mesh _grid, Parameters _param, std::vector<double>& conditions)
:m_grille(_grid), m_param(_param){

double dt = _grid.getdt();
double dx = _grid.getdx();
std::vector<double> m_init_vector = _grid.get_init_vector();
double v = _param.Get_Vol();
double r = _param.Get_Rate();
double theta = _param.Get_Theta();
long nb_step_time = _grid.Getvector_time().size();
long nb_step_spot = _grid.Getvector_stock().size();

//create the matrix of conditions with upper and lower filled and 0 otherwise
std::vector<double> cond(conditions);
std::vector<double> mat_conditions = transform_matrix(cond, nb_step_spot); 
        
//Sets the initial coefficients to define the matrices of the system to solve
double gamma_coefficient_previous = dt*theta*((v*v)/(dx*dx) + r); //creates the Gamma coefficient of index n+1
double gamma_coefficient = -gamma_coefficient_previous + dt; //creates the Gamma coefficient of index n
double alpha_coefficient_previous =  dt*theta*((-v*v)/(2*dx*dx) + (v*v)/(4*dx*dx) + (-r)/(2*dx)); //creates the Alpha coefficient of index n+1
double alpha_coefficient = -alpha_coefficient_previous + dt; // creates the Alpha coefficient of index n
double beta_coefficient_previous  =  dt*theta*((-v*v)/(2*dx*dx) + (v*v)/(4*dx*dx) + (r)/(2*dx)); //creates the Beta coefficient of index n+1
double beta_coefficient = -beta_coefficient_previous + dt; //creates the Beta coefficient of index n
           
//Constructs the matrices of dimension N-1 containing ones respectively on the lower diagonal, the diagonal, and the upper diagonal
xt::xarray<double> lower_diago = xt::eye(nb_step_spot,-1);
xt::xarray<double> diago = xt::eye(nb_step_spot);
xt::xarray<double> upper_diago = xt::eye(nb_step_spot,1);
            
//Constructs the initial values of the matrices A and B of the system to solve, thanks to the previously created matrices
xt::xarray<double> A = diago*(gamma_coefficient+1)+ lower_diago*(beta_coefficient) + upper_diago*(alpha_coefficient);
xt::xarray<double> B = diago*(-gamma_coefficient_previous+1) + lower_diago*(-beta_coefficient_previous) + upper_diago*(-alpha_coefficient_previous);

//Defines the vector B.X of the system to solve
xt::xarray<double> B_X;
    
//Converts the matrix of boundaries conditions into an xt::array type
xt::xarray<double> C_n = xt::adapt(cond);
        
//Converts the vector X_temp into an xt::array type vector named X
xt::xarray<double> X = xt::adapt(m_init_vector);
    

for (int t = 2; t <= nb_step_time; ++t){
    
        //Updates the time dependant coefficients
        gamma_coefficient_previous = gamma_coefficient;
        alpha_coefficient_previous = alpha_coefficient;
        beta_coefficient_previous = beta_coefficient;
        gamma_coefficient = -gamma_coefficient_previous + dt;
        alpha_coefficient = -alpha_coefficient_previous + dt;
        beta_coefficient = -beta_coefficient_previous + dt;
        
        //Updates the matrices A, B and B_X
        A = diago*(gamma_coefficient+1) + lower_diago*(beta_coefficient) + upper_diago*(alpha_coefficient);
        B = diago*(-gamma_coefficient_previous+1) + lower_diago*(-beta_coefficient_previous) + upper_diago*(-alpha_coefficient_previous);
		B_X = xt::linalg::dot(B,X);
    
        //Forces B_X as a vector
        B_X.reshape({nb_step_spot,static_cast<long>(1)});
    
        //Determines the value of X_(n) by solving the system: AX_(n) = B.X_(n+1) + C_(n+1) - C_(n)
        X = xt::linalg::solve(A, B_X + xt::view(C_n, xt::all(), xt::range(nb_step_time-t+1,nb_step_time-t+2)) - xt::view(C_n, xt::all(), xt::range(nb_step_time-t,nb_step_time-t+1)));


};

xt::xarray<double> _FX_n = X;

};

xt::xarray<double> Solve::Get_FX_n() const {
	return _FX_n;
};

Solve::~Solve(){}; //destructor of the solve object 


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Poubelle

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// PDE Solver

	// PDE::PDE(PayOff* _option, const double& dx,const double& dt,const std::vector<double>& time_vector,const std::vector<double>& spot_vector)
	// :option(_option)
	// {
		// size_t nb_step = spot_vector.size();
		// //std::vector<double> m_init_vector(nb_step);
		// for (std::size_t i = 0; i < nb_step ; ++i)
		// {
			// m_init_vector.push_back(init_cond(exp(spot_vector[i])));
		// }
		
	// };

	// // Initial condition (vanilla call option), we compute just the payoff created 
	// // x parameter stands for the spot
	// double PDE::init_cond(const double& x) const 
	// {
	  // return option->operator()(x);
	// };
	// std::vector<double> PDE::get_init_vector() const //const forbid to modify the state of my object
    // {
		// //m_nb_rows = 0;
        // return m_init_vector;
    // };


	// //Destructor
	// PDE::~PDE() {};


	

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
		