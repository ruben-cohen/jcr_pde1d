

#ifndef Resolve_HPP
#define Resolve_HPP

#include "mesh_spot.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace project{	
	
	class solver
	{
		public:
		
		solver(mesh grid, double theta,const std::vector<double>& boundaries, std::vector<std::vector<double>> vol_mat, std::vector<std::vector<double>> rate_mat); //constructor of the solver object
		
		~solver(); //destructor of the solver object 
		
		std::vector<double> Mid_diag_coeff(mesh grid,bool A,double theta,std::vector<double> sigma, std::vector<double> rate); 
		std::vector<double> Upper_diag_coeff(mesh grid,bool A,double theta,std::vector<double> sigma, std::vector<double> rate);
		std::vector<double> Lower_diag_coeff(mesh grid,bool A,double theta,std::vector<double> sigma, std::vector<double> rate);


		void thomas_algorithm(const std::vector<double>& upper_diag, const std::vector<double>& mid_diag, const std::vector<double>& lower_diag, const std::vector<double>& f_n1,
		std::vector<double>& f_sol);
		
		std::vector<std::vector<double>> get_price();
		
		std::vector<double> BX_vector(std::vector<double> upper, std::vector<double> mid, std::vector<double> low, std::vector<double> bound_diff,std::vector<double> Fn1);
		//this function will be used to compute at each time step the BX vector 
		
		private:
		
		mesh m_mesh;
		std::vector<std::vector<double>> m_results;
		
		//set of function to define the A matrix (3 better than just one huge matrix ?) 
		//in each function we only need to input the grid and the parameters as we will get theta, sigma and r from parameters class 
		// we will get dx and dt from grid class 
	};
}
#endif