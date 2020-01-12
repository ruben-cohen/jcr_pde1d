#include <iostream>
#include "closed_form.hpp"

int main(int argc, char* argv[])
{
	// Create the option parameters
	double K = 100;  // Strike price
	double r = 0.05;   // Risk-free rate (5%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 1.00;    // One year until expiry
	double res = 0.;
	// FDM discretisation parameters
	//double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
	//unsigned long J = 20; 
	//double t_dom = T;         // Time period as for the option
	//unsigned long N = 20;     

	// Create the PayOff and Option objects
	PayOff* pay_off_call = new PayOffCall(K);
	VanillaOption* call_option = new VanillaOption(K, r, T, v, pay_off_call);

	// Create the PDE and FDM objects
	//BlackScholesPDE* bs_pde = new BlackScholesPDE(call_option);
	//FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, bs_pde);
	//std::cout << *pay_off_call << std::endl;
	std::cout << "Julien" << std::endl;
	res = vanilla_payoff(100,90,1);
	std::cout << res << std::endl;
	//delete bs_pde;
	delete call_option;
	delete pay_off_call;

	return 0;
	
}
