#include 'Resolve.hpp'
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xadapt.hpp"
#include <vector>
#include <cmath>
using namespace xt::placeholders;

std::vector<double> resolution(mesh_spot grid, parameters param){
    
    //set the necessary parameters to define the matrices of the system to solve
    //double dt = grid.getdt(); //to get the time step
    //double dx = grid.getdx(); //to get the stock value step
    //double nb_time_steps = grid.Getvector_time.size(); //to get the number of time steps on which we iterate
    //double sigma = param.get_Vol(); //to get the volatility
    //double rate = param.get_Rate(); //to get the rate
    //double theta = param.get_Theta(); //to get the theta
    double size_matrix = nb_step - 1; //to get the size of the matrix A and B (N-1)
    
    //set the initial coefficients to define the matrices of the system to solve
    double gamma_coefficient_previous = dt*theta*((sigma**2)/(dx**2) + rate); //creates the Gamma coefficient of index n+1
    double gamma_coefficient_ = -gamma_coefficient_previous + dt; //creates the Gamma coefficient of index n
    double alpha_coefficient_previous =  dt*theta*((-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx)); //creates the Alpha coefficient of index n+1
    double alpha_coefficient = -alpha_coefficient_previous + dt; // creates the Alpha coefficient of index n
    double beta_coefficient_previous  =  dt*theta*((-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx)); //creates the Beta coefficient of index n+1
    double beta_coefficient = -beta_coefficient_previous + dt; //creates the Beta coefficient of index n
    
    //construct the matrices of dimension N-1 containing ones respectively on the lower diagonal, the diagonal, and the upper diagonal
    xt::xarray<double> lower_diago = xt::eye(size_matrix,-1);
    xt::xarray<double> diago = xt::eye(size_matrix);
    xt::xarray<double> upper_diago = xt::eye(size_matrix,1);
    
    //construct the initial values of the matrices A and B of the system to solve, thanks to the previously created matrices
    xt::xarray<double> A = diago*(gamma_coefficient+1) + lower_diago*(beta_coefficient) + upper_diago*(alpha_coefficient);
    xt::xarray<double> B = diago*(-gamma_coefficient_previous+1) + lower_diago*(-beta_coefficient_previous) + upper_diago*(-alpha_coefficient_previous);

    //construct the vectors C_n and C_n+1 of boundaries conditions --> Should be linked to the boundaries conditions class
    //std::vector<double> C_n(size_matrix); //Value taken from boundaries conditions
    //std::vector<double> C_n_previous(size_matrix); //Value taken from boundaries conditions
    
    //converts the matrix of boundaries conditions into an xt::array type
    xt::xarray<double> C_n = xt::adapt(conditions);
    
    //construct the vector X_temp to be determined by solving the system of equations. We initialize it to the Final payoff of the option
    std::vector<double> X_temp = m_init_vector;
    
    //converts the vector X_temp into an xt::array type vector named X
    xt::xarray<double> X = xt::adapt(X_temp);
    
    // We iterate on each time step of the grid to solve the equation. Beginning from time T-1 and solving for time 1.
    for (std::size_t t = 2; t < nb_time_steps; ++t){
    
        //we start by updating the time dependant coefficients
        gamma_coefficient_previous = gamma_coefficient;
        alpha_coefficient_previous = alpha_coefficient;
        beta_coefficient_previous = beta_coefficient;
        gamma_coefficient = -gamma_coefficient_previous + dt;
        alpha_coefficient = -alpha_coefficient_previous + dt;
        beta_coefficient = -beta_coefficient_previous + dt;
        
        //we then update the matrices A and B
        A = diago*(gamma_coefficient+1) + lower_diago*(beta_coefficient) + upper_diago*(alpha_coefficient);
        B = diago*(-gamma_coefficient_previous+1) + lower_diago*(-beta_coefficient_previous) + upper_diago*(-alpha_coefficient_previous);
    
        //We finally determine the value of X_(n) by solving the system: AX_(n) = BX_(n+1) + C_(n+1) - C_(n)
        X = xt::linalg::solve(A, xt::linalg::dot(B, X) + xt::view(C_n, xt::all(), xt::range(nb_time_steps-t-1,nb_time_steps-t)) - xt::view(C_n, xt::all(), xt::range(nb_time_steps-t-2,nb_time_steps-t-1)));

    };
    
    
};       
