


#ifndef BOUND_CONDITIONS_HPP
#define 
#include <vector>

// we create a boundaries conditions class that is purely virtual and 2 other classes neumann and derichtlet


namespace project{
	
	class bound_conditions {
	public: // all virtuals ? and we have all through neumann et/ou derichtlet ? 
	
	bound_conditions();
	virtual std::vector<double> operator() bounds(mesh_spot grid, parameters param, Payoff *function_payoff, std::vector<double> K_neuman ={0,0,0,0}) = 0;
	//the vector containing the neuman coef is passed as an optional parameters to avoid creating a new overload for neuman. Pas certain
	// que ce soit la meilleure solution on peut peut-être faire différement ? 
	//to be updated with payoff class from Julien 
	~bound_conditions();
	
	const std::vector<double> get_right_border();
	const std::vector<double> get_left_border();
	
	std::vector<double>  boundaries_compute(mesh_spot grid, parameters param, Payoff *function_payoff, bound_conditions*bound_func, std::vector<double> K_neuman ={0,0,0,0});
	//this function takes the same parameters as the bound_conditions + the type of boundaries we want to compute and will return the appropriate boundaries.

	private:
	
	std::vector<double> right_conditions;
	std::vector<double> left_conditions;
		
	};
	
	class Neumann public bound_conditions {
		
	neumann_conditions();
		
	virtual std::vector<double> operator() neumann(mesh_spot grid, parameters param, Payoff *function_payoff,std::vector<double> K_neuman ={0,0,0,0});
		
	
	}
	
	class Derichtlet public bound_conditions {
	
	Derichtlet_conditions();
	virtual std::vector<double> operator() derichtlet(mesh_spot grid, parameters param, Payoff *function_payoff,std::vector<double> K_neuman ={0,0,0,0});
		
		
		
		
		
	}
		
}


#endif 