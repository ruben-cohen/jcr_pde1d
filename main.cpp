#include <iostream>
#include "closed_form.hpp"
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
	
	// mesh discretisation parameters
	size_t nb_step_spot =20;    // Spot goes from [0.0, 1.0]
	long nb_step_time = 20; 
	
	// Create the PayOff object (strike as parameter)
	PayOff* option = new PayOffCall(K);
	// Create the mesh object (parameters (S,Vol, Nb time step, nb stock step)
	mesh grille(S,T,v,nb_step_time,nb_step_spot);
	//dauphine::matrix m(2, 4);

	
	//Create PDE object (to solve the problem, we provide to the class, the grid previoulsy defined in the class PDE) 
	PDE m(option,grille.getdx(),grille.getdt(),grille.Getvector_time(),grille.Getvector_stock());
	std::cout << m.init_cond(S) << std::endl;
	//print(m.get_init_vector());
	//print(grille.Getvector_time());
	m.print(m.get_init_vector());
	std::cout << grille.getdx() << std::endl;
	//delete bs_pde;
	delete option;

	return 0;
	
}
