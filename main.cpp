#include <iostream>
#include "closed_form.hpp"
#include <vector>

////
int main(int argc, char* argv[])
{
	// Create the option parameters
	double S = 100;
	double K = 100;  // Strike price
	double r = 0.05;   // Risk-free rate (5%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 1.00;    // One year until expiry
	double theta_ = 0.5;
	
	// mesh discretisation parameters
	long nb_step_spot =20;    // Spot goes from [0.0, 1.0]
	long nb_step_time = 20; 
	
	// Create the PayOff object (strike as parameter)
	PayOff* option = new PayOffCall(K);
	// Create the mesh object (parameters (S,Vol, Nb time step, nb stock step)
	mesh grille(S,T,v,nb_step_time,nb_step_spot);
	
	//Create PDE object (to solve the problem, we provide to the class, the grid previoulsy defined in the class PDE) 
	PDE m(option,grille.getdx(),grille.getdt(),grille.Getvector_time(),grille.Getvector_stock());
	
	//test of the boundaries : 
	Parameters par(v, r, theta_);
	std::vector<double> K_v{0.1,0.1,0.1,0.1};

	Derichtlet c(m,grille,par,option);
	Neumann c2(m,grille,par,option, K_v);
	
	std::cout<< "fonction calcul cond init:" << std::endl;
	std::cout << m.init_cond(S) << std::endl;
	std::cout<< "Vecteur de prix (log):" << std::endl;
	print(grille.Getvector_stock());
	std::cout<< "Vecteur de temps:" << std::endl;
	print(grille.Getvector_time());
	std::cout<< "Vecteur cond init (log):" << std::endl;
	print(m.get_init_vector());
	std::cout<< "dx:" << std::endl;
	std::cout << grille.getdx() << std::endl;
	std::cout<< "dt:" << std::endl;
	std::cout << grille.getdt() << std::endl;
	std::cout<< "vecteur dirichlet:" << std::endl;
	print(c.get_cond());
	std::cout<< "vecteur Neumann:" << std::endl;
	print(c2.get_cond());
	//delete bs_pde;
	delete option;
	// delete neu;
	//delete der;

	return 0;
	
}
