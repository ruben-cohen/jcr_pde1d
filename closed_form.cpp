
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

	std::vector<std::vector<double>> transform_matrix(std::vector<double> vector_init, double nb_rows){
	
	double endv; 
	//std::vector<double> upper_bound(vector_init);
	//std::vector<double> lower_bound(vector_init);
	
	
	 
	 if (vector_init.size()%2==0){endv=vector_init.size()/2;}
	 else {endv=std::floor(vector_init.size()/2) -1;}; //case useless !
	 
	 std::vector<double> upper_bound(endv);
	 std::vector<double> lower_bound(endv);
	 
	 std::copy(vector_init.begin(), vector_init.begin() + endv , upper_bound.begin());
	 
	 std::copy(vector_init.begin() + endv, vector_init.end(), lower_bound.begin());
	 
	 
	std::vector<double> row_0(upper_bound.size(), 0.0);
	 
	std::vector<std::vector<double>> matrix;
	
	matrix.resize(nb_rows, std::vector<double>(upper_bound.size()));
	
	matrix.front() = upper_bound;
	
	//std::cout << "row 0 size " << row_0.size() <<  "up size " << upper_bound.size() <<  "low size " << lower_bound.size() << std::endl;
	
	for(int i=1; i<nb_rows-1; i++){
		
		
		matrix[i] = row_0;
	};
	
	matrix.back() = lower_bound;
	 
	return matrix;		
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//parameters 
/* 
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

	}; */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Boundaries	
	
	bound_conditions::bound_conditions(mesh _grid)
	:
	 m_grille(_grid)
	 //,m_param(_param)
	 {};
	 
	
	Derichtlet::Derichtlet(const mesh& m_grid,const std::vector<double>& rate)
	:
	 bound_conditions(m_grid)
	{
		//From mesh object
		double dt = m_grid.getdt();
		double dx = m_grid.getdx(); //need the stock step 
		size_t size_vec = m_grid.Getvector_time().size();
		double maturity = m_grid.Getvector_time().back();
		double S0 = m_grid.get_Spot();
		size_t size_spot = m_grid.Getvector_stock().size();
		
/* 		if(size_spot != rate.size()){
			
			std::cout<< "size of sigma is" << rate.size() << "size of spot is" << size_spot << std::endl;
		}
		else{ std::cout<< "okay for Derichtlet" << std::endl;}; */
		
		//From parametre object
		//double sigma = m_param.Get_Vol(); //to get the volatility 
		//double rate = m_param.Get_Rate(); //the get the rate 
		//double theta = m_param.Get_Theta(); //to get the theta 
		//double r = m_param.Get_Rate();
		
		//From PDE object
		std::vector<double>  _init_cond = m_grid.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		double f_0_T = _init_cond[0];
		double f_N_T = _init_cond.back();
		
		matrix_derichtlet.resize(size_vec*2); 
		//std::vector<double> lower_conditions(size_vec); 
		
		for(size_t i = 0; i < size_vec; ++i)
        {
            matrix_derichtlet[i] = f_0_T*exp(-rate[i]*dt*i); // for here
			matrix_derichtlet[size_vec+i] = f_N_T*exp(-rate[i]*dt*i);
        }
	};
	
	std::vector<double> Derichtlet::get_cond() const
	{
		return matrix_derichtlet;
	};
	
	Neumann::Neumann(mesh m_grid, double theta, std::vector<double> sigma, std::vector<double> rate,std::vector<double>& const_vector)
	:bound_conditions(m_grille)
	{
		//From mesh object
		double dt = m_grid.getdt();
		double dx = m_grid.getdx(); //need the stock step 
		size_t size_vec = m_grid.Getvector_time().size();
		double maturity = m_grid.Getvector_time().back();
		double S0 = m_grid.get_Spot();
		size_t size_spot = m_grid.Getvector_stock().size();		
		
		//From parametre object
		//double sigma = m_param.Get_Vol(); //to get the volatility 
		//double rate = m_param.Get_Rate(); //the get the rate 
		//double theta = m_param.Get_Theta(); //to get the theta 
		//double r = m_param.Get_Rate();
		
/* 		if(size_spot != sigma.size()){
			
			std::cout<< "size of sigma is" << sigma.size() << "size of spot is" << size_spot << std::endl;
		}
		else{ std::cout<< "okay for neumann" << std::endl;};
		
		if(size_spot != rate.size()){
			
			std::cout<< "size of rate is" << rate.size() << "size of spot is" << size_spot << std::endl;
		}
		else{ std::cout<< "okay for neumann" << std::endl;}; */
		
		//From PDE object
		std::vector<double>  _init_cond = m_grid.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		double f_0_T = _init_cond[0];
		double f_N_T = _init_cond.back();
		
	    double K1 = const_vector[0];
	    double K2 = const_vector[1];
	    double K3 = const_vector[2];
	    double K4 = const_vector[3];
		   
		double coef_; 
		double coef_K1_K2; 
		double coef_K3_K4; 
		
		matrix_neumann.resize(size_vec);
		matrix_neumann.back() =  f_0_T*exp(-maturity*rate[0]);
		matrix_neumann.resize(size_vec*2);
		matrix_neumann.back() = f_N_T*exp(-maturity*rate.back());
		
		for (size_t it = size_vec-1; it > 0; it--)
		{
		      //reverse iterator to fill the vector from the end to the beginning
			  
			  //std::cout << it << std::endl;
			  
			  size_t tt = size_vec + it -1;
			  
			  //std::cout << tt << std::endl;
			  
			  		coef_ =(1 - dt*(1-theta)*rate[it])/(1+ dt*theta*rate[it]);
					coef_K1_K2 = -dt*((-std::pow(sigma[it],2)*K1)/2 + (std::pow(sigma[it],2)/2 - rate[it])*K2)/(1+ dt*theta*rate[it]);
					coef_K3_K4 = -dt*((-std::pow(sigma[it],2)*K3)/2 + (std::pow(sigma[it],2)/2 - rate[it])*K4)/(1+ dt*theta*rate[it]);
					
					
			  matrix_neumann[it]  = (coef_*matrix_neumann[it-1] + coef_K1_K2)*exp(-rate[it]*dt*it);
			
			  matrix_neumann[tt] = (coef_*matrix_neumann[tt-1] + coef_K3_K4)*exp(-rate[it]*dt*it);
		 }
		 
	};
	
    std::vector<double> Neumann::get_cond() const
	{
		return matrix_neumann;
	};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Classe de base
	//Constructor
	volatility::volatility(const double& v,const mesh& grid)
	:
	 m_init_vol(v),
	 m_grid(grid)
	{
		m_vector_time = m_grid.Getvector_time();
		m_vector_stock = m_grid.Getvector_stock();
		m_nb_time = m_vector_time.size();
		m_nb_spot = m_vector_stock.size();	
	};
	
	//Here no computation since vol is constant --> just return init vol 
	double volatility::compute_vol(const double& S,const double& t,const double& x,const double& y)
	{
		return m_init_vol;
	};
	
	std::vector<double> volatility::vector_vol(const std::vector<double> v_spot,const double& t)
	{
		std::size_t nb_step = v_spot.size();
		
		//m_vol_const.resize(nb_step,0.);
		
		for (long i = 0; i < nb_step ; ++i)
		{
			m_vol_const.push_back(compute_vol());
		}
		
		return m_vol_const;
	};
//////////////////////////////////////////////////////////////////////////////
//Classe qui hérite pour la surface de vol

	vol_surface::vol_surface(const double& v, mesh grid,const double& coeff_tps, const double& coeff_spot)
	:
	 volatility(v, grid), //Normalament on obtient ainsi les variables protected de volatility
	 m_coeff_spot(coeff_spot),
	 m_coeff_time(coeff_tps)
	{
	};
	
	//Here computation of a vol thanks to parameters
	double vol_surface::compute_vol(const double& S,const double& t,const double& x,const double& y)
	{
		double vol = x*t + y*S;
		
		return vol;
	};
	
	std::vector<double> vol_surface::vector_vol(const std::vector<double> v_spot,const double& t)
	{
		std::size_t nb_step = v_spot.size();
		
		//m_vol_const.resize(nb_step,0.);
		
		for (long i = 0; i < nb_step ; ++i)
		{
			m_vol_matrix.push_back(compute_vol(v_spot[i],t,m_coeff_time,m_coeff_spot));
		}
		
		return m_vol_matrix;
	};	



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 solver::solver(mesh grid, std::vector<std::vector<double>> res):m_mesh(grid),m_results(res){};
 
 void solver::solve_X(mesh grid, double theta, const std::vector<double>& boundaries, std::vector<std::vector<double>> vol_mat,std::vector<std::vector<double>> rate_mat)
 {
	 
	 
	 std::vector<std::vector<double>> results;
	 std::vector<std::vector<double>> vol(vol_mat);
	 std::vector<std::vector<double>> rate(rate_mat);
	 long T = grid.Getvector_time().size(); 
	 std::size_t m = grid.Getvector_stock().size();
	 std::vector<double> init_f = grid.get_init_vector(); // get the initial X_T
	 std::vector<std::vector<double>> matrix_boundaries = transform_matrix(boundaries,T); //get the boundaries matrix 
	 std::vector<double> cond_diff(m);
	 std::vector<double> solution_vector(m-2, 0.0); //init of the solution vector 
	 std::vector<double> up_A;
	 std::vector<double> low_A;
	 std::vector<double> mid_A;
	 std::vector<double> up_B;
	 std::vector<double> low_B;
	 std::vector<double> mid_B;
	 
	 solver solver(grid, results);
	 
	 
	 for(size_t s = 0; s<m; ++s)
	 {
		 
		 cond_diff[s] = matrix_boundaries[s][T]-matrix_boundaries[s][T-1];
	 }; //initialisation du vecteur de C(n+1) - C(n) 
	 
	 //initialisation de la matrice B à maturité 
	 up_B = solver.Upper_diag_coeff(grid,false,theta,vol.back(), rate.back());
	 low_B = solver.Lower_diag_coeff(grid,false,theta,vol.back(), rate.back());
	 mid_B = solver.Mid_diag_coeff(grid,false,theta,vol.back(), rate.back());
	 
	 print(up_B);
	 print(low_B);
	 print(mid_B);
	 print(init_f);
	 
	 
	std::vector<double> B = solver.BX_vector(up_B,mid_B,low_B,cond_diff,init_f);

	rate.pop_back();
	vol.pop_back();
	 
	 
	 for(int i = grid.Getvector_time().size() - 1; i > 0; --i)
	 {
		 
		
		up_A = solver.Upper_diag_coeff(grid,true,theta,vol.back(), rate.back());
		low_A = solver.Lower_diag_coeff(grid,true,theta,vol.back(), rate.back());
		mid_A = solver.Mid_diag_coeff(grid,true,theta,vol.back(), rate.back());
		
		print(up_A);
		print(low_A);
		print(mid_A);
			
		solver.thomas_algorithm(up_A, mid_A, low_A, B, solution_vector);
		
		print(solution_vector);
		results.push_back(solution_vector); //create a matrix containing all the prices computed by the solver in order to be displayed afterwards
		
		for(size_t s = 0; s<m; ++s){
		 
				cond_diff[s] = matrix_boundaries[s][i]-matrix_boundaries[s][i-1];
			 }; 
			 
		print(cond_diff);
			 
		//init_f = solution_vector;
		
		up_B = solver.Upper_diag_coeff(grid, false,theta,vol.back(), rate.back());
	    low_B = solver.Lower_diag_coeff(grid, false,theta,vol.back(), rate.back());
		mid_B = solver.Mid_diag_coeff(grid, false,theta,vol.back(), rate.back());
	 
	 
	    B = solver.BX_vector(up_B,mid_B,low_B,cond_diff,solution_vector);
		
		rate.pop_back();
		vol.pop_back();
		
		 
	 }; 
	  
 }; 
 
std::vector<std::vector<double>> solver::get_price(){
	
	return m_results;
};
 
 std::vector<double> solver::BX_vector(std::vector<double> upper, std::vector<double> mid, std::vector<double> low, std::vector<double> bound_diff, std::vector<double> Fn1){
//this procedure creates the right_hand vector of the AX = D equation based on BX(n+1) + C(n+1) - C(n)
	std::size_t N = mid.size(); // as the resolution si between f1 and fN-1
	
	std::vector<double> BX(N);
	
	for(size_t i = 0; i < N; i++){
		
		
		if(i==0){BX[0] = mid[0]*Fn1[0] + upper[0]*Fn1[0] +bound_diff[0];}
			
			
		else if(i==N) {BX[N] =  mid.back()*Fn1.back() + low.back()*Fn1[N-1] + bound_diff.back();}
        
    
	
		else{BX[i] = mid[i]*Fn1[i] + upper[i]*Fn1[i+1] +low[i-1]*Fn1[i-1] + bound_diff[i];}
	
	};
	
	return BX;
};
	
//set of function to define the A matrix (3 better than just one huge matrix ?) 
std::vector<double> solver::Mid_diag_coeff(mesh grid, bool A,double theta, std::vector<double> sigma, std::vector<double> rate){
	
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	//std::vector<double> sigma_init = param.get_Vol(); //to get the volatility 
	//std::vector<double> rate_init = param.get_Rate(); //the get the rate 
	//double theta = param.get_Theta(); //to get the theta 
	double size_gamma = grid.Getvector_stock().size(); //no minus-1 as we are on diagonal 
	
	
	
		//create the vector that holds the diagonal 
	std::vector<double> gamma_coefficient(size_gamma);
	
	for (std::size_t i = 0; i < size_gamma; ++i){
		
		gamma_coefficient[i] =  std::pow(sigma[i],2)/std::pow(dx,2) + rate[i];
		
		if (A==false){
			
			gamma_coefficient[i] = -dt*(1-theta)*gamma_coefficient[i] + 1;
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		}
		
		else{
			
			gamma_coefficient[i] = dt*theta*gamma_coefficient[i] + 1;
		};
	};
	
	return gamma_coefficient;
	
}; 


std::vector<double> solver::Upper_diag_coeff(mesh grid, bool A,double theta,std::vector<double> sigma, std::vector<double> rate){
	
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	//std::vector<double> sigma_init = param.get_Vol(); //to get the volatility 
	//std::vector<double> rate_init = param.get_Rate(); //the get the rate 
	//double theta = param.get_Theta(); //to get the theta 
	double size_alpha = grid.Getvector_stock().size()-1; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> alpha_coefficient(size_alpha);
	
	for (std::size_t i = 0; i < size_alpha; ++i){
		
		alpha_coefficient[i] =  (-std::pow(sigma[i],2))/(2*std::pow(dx,2)) + std::pow(sigma[i],2)/(4*std::pow(dx,2)) + (-rate[i])/(2*dx);
		
		 
		if (A==false){
			
			alpha_coefficient[i] = -dt*(1-theta)*alpha_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		}
		
		else{
			
			alpha_coefficient[i] = dt*theta*alpha_coefficient[i];
		};
	};
	
	return alpha_coefficient;
	
};


std::vector<double> solver::Lower_diag_coeff(mesh grid, bool A,double theta,std::vector<double> sigma, std::vector<double> rate){
	
	
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	//std::vector<double> sigma_init = param.get_Vol(); //to get the volatility 
	//std::vector<double> rate_init = param.get_Rate(); //the get the rate 
	//double theta = param.get_Theta(); //to get the theta 
	double size_beta = grid.Getvector_stock().size()-1; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> beta_coefficient(size_beta);
	
	for (std::size_t i = 0; i < size_beta; ++i){
		
		beta_coefficient[i] =  (-std::pow(sigma[i],2))/(2*std::pow(dx,2)) + std::pow(sigma[i],2)/(4*std::pow(dx,2)) + (rate[i])/(2*dx);
		
		
		if (A==false){
			
			beta_coefficient[i] = -dt*(1-theta)*beta_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		}
		
		else{
			
			beta_coefficient[i] = dt*theta*beta_coefficient[i];
		};
	};
	
	return beta_coefficient;
	
};

//this is the thomas algo for inverting the matrix  

void solver::thomas_algorithm(const std::vector<double>& upper_diag,
                      const std::vector<double>& mid_diag,
                      const std::vector<double>& lower_diag,
                      const std::vector<double>& f_n1,
                      std::vector<double>& f_sol) {
  size_t nb_spot = f_n1.size();

  // Create the temprary vectors to store new coef                                                                                                                                                                                                                                                                                                                                                
  std::vector<double> c_star(nb_spot, 0.0);
  std::vector<double> d_star(nb_spot, 0.0);                                                                                                                                                    
  c_star[0] = lower_diag[0] / mid_diag[0];
  d_star[0] = f_n1[0] / mid_diag[0];

  //forward sweep                                                                                                                                                  
  for (int i=1; i<nb_spot; i++) {
    double m = 1.0 / (mid_diag[i] - upper_diag[i] * c_star[i-1]);
    c_star[i] = lower_diag[i] * m;
    d_star[i] = (f_n1[i] - upper_diag[i] * d_star[i-1]) * m;
  }

  //reverse sweep, used to update the solution vector f                                                                                                                                                 
  for (int i=nb_spot-1; i-- > 0; ) {
    f_sol[i] = d_star[i] - c_star[i] * f_n1[i+1];
  }
}; 


solver::~solver(){};

//this is a void function because it modifies the price 
	
/*  Solve::Solve(mesh _grid, Parameters _param, std::vector<double>& conditions)
 :
	m_grid(_grid),
	m_param(_param)
 {

	 double dt = m_grid.getdt();
	 double dx = m_grid.getdx();
	 std::vector<double> m_init_vector = m_grid.get_init_vector();
	 double v = m_param.Get_Vol();
	 double r = m_param.Get_Rate();
	 double theta = m_param.Get_Theta();
	 long nb_step_time = m_grid.Getvector_time().size();
	 long nb_step_spot = m_grid.Getvector_stock().size();

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
		

	 for (int t = 2; t <= nb_step_time; ++t)
	{
		
			// Updates the time dependant coefficients
			gamma_coefficient_previous = gamma_coefficient;
			alpha_coefficient_previous = alpha_coefficient;
			beta_coefficient_previous = beta_coefficient;
			gamma_coefficient = -gamma_coefficient_previous + dt;
			alpha_coefficient = -alpha_coefficient_previous + dt;
			beta_coefficient = -beta_coefficient_previous + dt;
			
			// Updates the matrices A, B and B_X
			A = diago*(gamma_coefficient+1) + lower_diago*(beta_coefficient) + upper_diago*(alpha_coefficient);
			B = diago*(-gamma_coefficient_previous+1) + lower_diago*(-beta_coefficient_previous) + upper_diago*(-alpha_coefficient_previous);
			B_X = xt::linalg::dot(B,X);
		
			// Forces B_X as a vector
			B_X.reshape({nb_step_spot,static_cast<long>(1)});
		
			// Determines the value of X_(n) by solving the system: AX_(n) = B.X_(n+1) + C_(n+1) - C_(n)
			X = xt::linalg::solve(A, B_X + xt::view(C_n, xt::all(), xt::range(nb_step_time-t+1,nb_step_time-t+2)) - xt::view(C_n, xt::all(), xt::range(nb_step_time-t,nb_step_time-t+1)));


	};

	xt::xarray<double> _FX_n = X;

 };

xt::xarray<double> Solve::Get_FX_n() const 
{
	//return _FX_n;
//};
 */
//Solve::~Solve(){}; //destructor of the solve object 


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
		