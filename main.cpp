#include <iostream>
#include "closed_form.hpp"

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
	
	// FDM discretisation parameters
	//double x_dom_ma = 1.0;       // Spot goes from [0.0, 1.0]
	//unsigned long J = 20; 
	//double t_dom = T;         // Time period as for the option
	//unsigned long N = 20;     
	//PayOff* pay_off_call = new PayOffCall(K);
	
// Create the PayOff object (strike as parameter)
	PayOff* option = new PayOffCall(K);
// Create the PDE objects as follows (the class should define the mesh given the FDM discretisation parameters)
	//PDE* bs_pde = new PDE(option);
	
//Create FDM object (to solve the problem, we provide to the class, the grid previoulsy defined in the class PDE) 
	
	//Test pour voir si la fonction de calcul initial condition fonctionne
	PDE m(option);
	std::cout << m.init_cond(S) << std::endl;
	
	//FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, bs_pde);

	
	std::cout << "Julien" << std::endl;
	//delete bs_pde;
	delete option;

	return 0;
	
}
