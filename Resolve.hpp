

#ifndef Resolve_HPP
#define Resolve_HPP

#include 'mesh_spot.hpp'
#include <vector>

namespace project {
	
	
	std::vector<double> CranckNicholson_algo(mesh_spot grid, parameters param); //main function where we define the CN procedure
	std::vector<double> Thomas_triLinMatrix_inverse(); //need to inverse A in the system that is a diagonal matrix 
	
	//set of function to define the A matrix (3 better than just one huge matrix ?) 
	//in each function we only need to input the grid and the parameters as we will get theta, sigma and r from parameters class 
	// we will get dx and dt from grid class 
	std::vector<double> Mid_diag_coeff(mesh_spot grid, parameters param,bool A); 
	std::vector<double> Upper_diag_coeff(mesh_spot grid, parameters param,bool A);
	std::vector<double> Lower_diag_coeff(mesh_spot grid, parameters param,bool A);
	
	//all in function to get the price of the option and the greeks ? 
    std::vector<double> resolution();
	
	//function to compute the greeks ? 
	
}






#endif