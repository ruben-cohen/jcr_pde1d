#ifndef Resolution_hpp
#define Resolution_hpp

#include <vector>
#include <cmath>
#include <iostream>
//Allows to create and manipulate xarrays
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xeval.hpp"
//Allows to perform linear algebra operations on xarrays
#include "xtensor-blas/xlinalg.hpp"

namespace project
{
    class Solve
    {
        
    public:

        Solve();
        
        double dt; //to get the time step
        double dx; //to get the stock value step
        int nb_step_time; //to get the number of time steps on which we iterate
        double v; //to get the volatility
        double r; //to get the rate
        double theta; //to get the theta
        int nb_step_spot; //to get the size of the matrix A and B (N-1)
        xt::xarray<double> C_n(std::vector<double> conditions); //to get the matrix of boundaries conditions
        xt::xarray<double> X(std::vector<double> m_init_vector); //to get the initial value of the vector X to solve
        
        //Sets the initial coefficients to define the matrices of the system to solve
        double gamma_coefficient_previous = dt*theta*((v*v)/(dx*dx) + r); //creates the Gamma coefficient of index n+1
        double gamma_coefficient = -gamma_coefficient_previous + dt; //creates the Gamma coefficient of index n
        double alpha_coefficient_previous =  dt*theta*((-v*v)/(2*dx*dx) + (v*v)/(4*dx*dx) + (-r)/(2*dx)); //creates the Alpha coefficient of index n+1
        double alpha_coefficient = -alpha_coefficient_previous + dt; // creates the Alpha coefficient of index n
        double beta_coefficient_previous  =  dt*theta*((-v*v)/(2*dx*dx) + (v*v)/(4*dx*dx) + (r)/(2*dx)); //creates the Beta coefficient of index n+1
        double beta_coefficient = -beta_coefficient_previous + dt; //creates the Beta coefficient of index n
           
        //Constructs the matrices of dimension N-1 containing ones respectively on the lower diagonal, the diagonal, and the upper diagonal
        xt::xarray<double> lower_diago = xt::eye(nb_step_spot,-1);
        xt::xarray<double> diago = xt::eye(nb_step_spot);
        xt::xarray<double> upper_diago = xt::eye(nb_step_spot,1);
            
        //Constructs the initial values of the matrices A and B of the system to solve, thanks to the previously created matrices
        xt::xarray<double> A = diago*(gamma_coefficient+1)+ lower_diago*(beta_coefficient) + upper_diago*(alpha_coefficient);
        xt::xarray<double> B = diago*(-gamma_coefficient_previous+1) + lower_diago*(-beta_coefficient_previous) + upper_diago*(-alpha_coefficient_previous);

        //Defines the vector B.X of the system to solve
        xt::xarray<double> B_X;

        //Iterates on each time step of the grid to solve the equation. Beginning from time T-1 and solving for time 1

    };
}
#endif /* Resolution_hpp */
