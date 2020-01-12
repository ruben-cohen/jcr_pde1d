
#include 'Resolve.hpp'
#include <vector>
#include <cmath>


PDE::PDE(VanillaOption* _option)
 : option(_option)
{
}

// Initial condition (vanilla call option)
double PDE::init_cond(double x) const 
{
  return option->pay_off->operator()(x);
}

std::vector<double> CranckNicholson_algo(mesh_spot grid, parameters param){}; //main function where we define the CN procedure
std::vector<double> Thomas_triLinMatrix_inverse(){}; //need to inverse A in the system that is a diagonal matrix 
	
//set of function to define the A matrix (3 better than just one huge matrix ?) 
std::vector<double> Mid_diag_coeff(mesh_spot grid, parameters param,bool A){
	
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	double sigma = param.get_Vol(); //to get the volatility 
	double rate = param.get_Rate(); //the get the rate 
	double theta = param.get_Theta(); //to get the theta 
	double size_gamma = Getvector_stock().size(); //no minus-1 as we are on diagonal 
	
		//create the vector that holds the diagonal 
	std::vector<double> gamma_coefficient(size_gamma);
	
	for (std::size_t i = 0; i < size_gamma; ++i){
		
		gamma_coefficient[i] =  (sigma**2)/(dx**2) + rate;
		
		if (A==false){
			
			gamma_coefficient[i] = -dt*(1-theta)*gamma_coefficient[i] + 1;
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		};
		
		else{
			
			gamma_coefficient[i] = dt*theta*alpha_coefficient[i] + 1;
		};
	};
	
	return gamma_coefficient;
	
}; 


std::vector<double> Upper_diag_coeff(mesh_spot grid, parameters param,bool A){
	
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	double sigma = param.get_Vol(); //to get the volatility 
	double rate = param.get_Rate(); //the get the rate 
	double theta = param.get_Theta(); //to get the theta 
	double size_alpha = Getvector_stock().size()-1; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> alpha_coefficient(size_alpha);
	
	for (std::size_t i = 0; i < size_alpha; ++i){
		
		alpha_coefficient[i] =  (-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
		
		 
		if (A==false){
			
			alpha_coefficient[i] = -dt*(1-theta)*alpha_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		};
		
		else{
			
			alpha_coefficient[i] = dt*theta*alpha_coefficient[i];
		};
	};
	
	return alpha_coefficient;
	
};


std::vector<double> Lower_diag_coeff(mesh_spot grid, parameters param,bool A){
	
	
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	double sigma = param.get_Vol(); //to get the volatility 
	double rate = param.get_Rate(); //the get the rate 
	double theta = param.get_Theta(); //to get the theta 
	double size_beta = Getvector_stock().size()-1; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> beta_coefficient(size_beta);
	
	for (std::size_t i = 0; i < size_beta; ++i){
		
		beta_coefficient[i] =  (-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
		
		
		if (A==false){
			
			beta_coefficient[i] = -dt*(1-theta)*beta_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		};
		
		else{
			
			beta_coefficient[i] = dt*theta*beta_coefficient[i];
		};
	};
	
	return beta_coefficient;
	
};
//all in function to get the price of the option and the greeks ? 
std::vector<double> resolution(){};
	
//function to compute the greeks ? 