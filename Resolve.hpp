

#ifndef Resolve_HPP
#define Resolve_HPP

#include <vector>

namespace project {
	
	
	std::vector<double> CranckNicholson_algo(); //main function where we define the CN procedure
	std::vector<double> Thomas_triLinMatrix_inverse(); //need to inverse A in the system that is a diagonal matrix 
	
	//set of function to define the A matrix (3 better than just one huge matrix ?) 
	std::vector<double> Mid_diag_coeff(double theta); 
	std::vector<double> Upper_diag_coeff(double theta);
	std::vector<double> Lower_diag_coeff(double theta);
	
	//all in function to get the price of the option and the greeks ? 
    std::vector<double> resolution();
	
	//function to compute the greeks ? 
	
	
	
	
	
	
	
	
	
	
	
	
	
}






#endif