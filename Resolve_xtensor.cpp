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



//ACTIVATE FOR TESTING ONLY
//Sets the necessary parameters to define the matrices of the system to solve --> arbitrary values used for tests
double dt = 0.1; //to get the time step
double dx = 0.1; //to get the stock value step
int nb_time_steps = 5; //to get the number of time steps on which we iterate
double sigma = 0.2; //to get the volatility
double rate = 0.01; //to get the rate
double theta = 0.5; //to get the theta
int size_matrix = 10; //to get the size of the matrix A and B (N-1)
xt::xarray<double> C_n = xt::eye(size_matrix,1); //to get the matrix of boundaries conditions
xt::xarray<double> X = xt::linspace<double>(1.0, size_matrix, size_matrix); //to get the initial value of the vector X to solve


//COMMENT THESE 3 LINES FOR TESTING ONLY
//Converts the matrix of boundaries conditions into an xt::array type
//xt::xarray<double> C_n = xt::adapt(conditions);
    
//Constructs the vector X_temp that we initialize to the final payoff of the option
//std::vector<double> X_temp = m_init_vector;

//Converts the vector X_temp into an xt::array type vector named X
//xt::xarray<double> X = xt::adapt(X_temp);



//Sets the initial coefficients to define the matrices of the system to solve
double gamma_coefficient_previous = dt*theta*((sigma*sigma)/(dx*dx) + rate); //creates the Gamma coefficient of index n+1
double gamma_coefficient = -gamma_coefficient_previous + dt; //creates the Gamma coefficient of index n
double alpha_coefficient_previous =  dt*theta*((-sigma*sigma)/(2*dx*dx) + (sigma*sigma)/(4*dx*dx) + (-rate)/(2*dx)); //creates the Alpha coefficient of index n+1
double alpha_coefficient = -alpha_coefficient_previous + dt; // creates the Alpha coefficient of index n
double beta_coefficient_previous  =  dt*theta*((-sigma*sigma)/(2*dx*dx) + (sigma*sigma)/(4*dx*dx) + (rate)/(2*dx)); //creates the Beta coefficient of index n+1
double beta_coefficient = -beta_coefficient_previous + dt; //creates the Beta coefficient of index n
   
//Constructs the matrices of dimension N-1 containing ones respectively on the lower diagonal, the diagonal, and the upper diagonal
xt::xarray<double> lower_diago = xt::eye(size_matrix,-1);
xt::xarray<double> diago = xt::eye(size_matrix);
xt::xarray<double> upper_diago = xt::eye(size_matrix,1);
    
//Constructs the initial values of the matrices A and B of the system to solve, thanks to the previously created matrices
xt::xarray<double> A = diago*(gamma_coefficient+1)+ lower_diago*(beta_coefficient) + upper_diago*(alpha_coefficient);
xt::xarray<double> B = diago*(-gamma_coefficient_previous+1) + lower_diago*(-beta_coefficient_previous) + upper_diago*(-alpha_coefficient_previous);

//Defines the vector B.X of the system to solve
xt::xarray<double> B_X;

//Iterates on each time step of the grid to solve the equation. Beginning from time T-1 and solving for time 1
for (int t = 2; t <= nb_time_steps; ++t){
    
        //Updates the time dependant coefficients
        gamma_coefficient_previous = gamma_coefficient;
        alpha_coefficient_previous = alpha_coefficient;
        beta_coefficient_previous = beta_coefficient;
        gamma_coefficient = -gamma_coefficient_previous + dt;
        alpha_coefficient = -alpha_coefficient_previous + dt;
        beta_coefficient = -beta_coefficient_previous + dt;
        
        //Updates the matrices A, B and B_X
        A = diago*(gamma_coefficient+1) + lower_diago*(beta_coefficient) + upper_diago*(alpha_coefficient);
        B = diago*(-gamma_coefficient_previous+1) + lower_diago*(-beta_coefficient_previous) + upper_diago*(-alpha_coefficient_previous);
        B_X = xt::linalg::dot(B, X);
    
        //Forces B_X as a vector
        B_X.reshape({size_matrix,1});
    
        //Determines the value of X_(n) by solving the system: AX_(n) = B.X_(n+1) + C_(n+1) - C_(n)
        X = xt::linalg::solve(A, B_X + xt::view(C_n, xt::all(), xt::range(nb_time_steps-t+1,nb_time_steps-t+2)) - xt::view(C_n, xt::all(), xt::range(nb_time_steps-t,nb_time_steps-t+1)));
};
