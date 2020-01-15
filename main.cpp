#include <iostream>
#include "closed_form.hpp"
#include <vector>

//
int main(int argc, char* argv[])
{
	// Create the option parameters
	double S = 105;
	double K = 100;  // Strike price
	double r = 0.05;   // Risk-free rate (5%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 1.00;    // One year until expiry
	double res = 0.;
	double res_2 = 0.;
	double theta_ = 0.5;
	
	// mesh discretisation parameters
	size_t nb_step_spot =20;    // Spot goes from [0.0, 1.0]
	long nb_step_time = 20; 
	
	// Create the PayOff object (strike as parameter)
	PayOff* option = new PayOffCall(K);
	// Create the mesh object (parameters (S,Vol, Nb time step, nb stock step)
	mesh grille(S,T,v,nb_step_time,nb_step_spot);
	//dauphine::matrix m(2, 4);
	
	//test of the boundaries : 
	Parameters par(v, r, theta_);
	std::vector<double> K_v{1.,2.,3.,4.};
	bound_conditions* neu = new Neumann();
	std::vector<double> K_d{0.,0.,0.,0.};
	bound_conditions* der = new Derichtlet();
	
	//Create PDE object (to solve the problem, we provide to the class, the grid previoulsy defined in the class PDE) 
	PDE m(option,grille.getdx(),grille.getdt(),grille.Getvector_time(),grille.Getvector_stock());
	
	std::vector<std::vector<double>> der_matrix  = bound_conditions::boundaries_compute(grille,par,option,der, K_d);
	
	std::vector<std::vector<double>> neu_matrix = bound_conditions::boundaries_compute(grille,par,option,neu, K_v);
	
	//std::cout< "Derichtlet conditions matrix" < std::endl;
	
	for(size_t i = 0; i< der_matrix.size(); i++){
		
		for (size_t j = 0; j < der_matrix[i].size(); j++){
		
		std::cout<< der_matrix[i][j] ;
		
		}
		
		std::cout<< " \n" << std::endl;
	}
	
	std::cout<< "Neumann conditions matrix" << std::endl;
	
	for(size_t i = 0; i<neu_matrix.size(); i++){
		
		for (size_t j = 0; j < neu_matrix[i].size(); j++){
		
		std::cout<<neu_matrix[i][j] ;
		
		}
		
		std::cout<< " \n" << std::endl;
	}
	
	std::cout << m.init_cond(S) << std::endl;
	//print(m.get_init_vector());
	//print(grille.Getvector_time());
	m.print(m.get_init_vector());
	std::cout << grille.getdx() << std::endl;
	//delete bs_pde;
	delete option;

	return 0;
	
}
