
#include 'Resolve.hpp'
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include <vector>
#include <cmath>


PDE::PDE(VanillaOption* _option)
 : option(_option)
{
}

// Initial condition (vanilla call option)
double PDE::init_cond(double x) const
{
  return option->pay_off->operator()(x);
}

std::vector<double> CranckNicholson_algo(mesh_spot grid, parameters param){}; //main function where we define the CN procedure


std::vector<double> matrices(mesh_spot grid, parameters param){
    
    //set the necessary parameters to define the matrices of the system to solve
    double dt = grid.getdt(); //need the time step
    double dx = grid.getdx(); //need the stock step
    double sigma = param.get_Vol(); //to get the volatility
    double rate = param.get_Rate(); //the get the rate
    double theta = param.get_Theta(); //to get the theta
    
    double size_matrix = Getvector_stock().size() - 1; //to get the size of the matrix A and B (N-1)
    double gamma_coefficient = (sigma**2)/(dx**2) + rate; //creates the Gamma variable of index n
    double gamma_coefficient_next = -gamma_coefficient + dt; //creates the Gamma variable of index n+1
    double alpha_coefficient =  (-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx); //creates the Alpha variable of index n
    double alpha_coefficient_next = -alpha_coefficient + dt; // creates the Alpha variable of index n+1
    double beta_coefficient  =  (-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx); //creates the Beta variable of index n
    double beta_coefficient_next = -beta_coefficient + dt; //creates the Beta variable of index n+1
    
    //construct the matrices of dimension N-1 containing ones respectively on the lower diagonal, the diagonal, and the upper diagonal
    xt::xarray<double> lower_diago = xt::eye(size_matrix,-1);
    xt::xarray<double> diago = xt::eye(size_matrix);
    xt::xarray<double> upper_diago = xt::eye(size_matrix,1);
    
    //construct the matrices A and B of the system to solve
    xt::xarray<double> A = diago*(gamma_coefficient+1) + lower_diago*(beta_coefficient) + upper_diago*(alpha_coefficient);
    xt::xarray<double> B = diago*(-gamma_coefficient_next+1) + lower_diago*(-beta_coefficient_next) + upper_diago*(-alpha_coefficient_next);

    //construct the vectors C_n and C_n+1 of boundaries conditions --> Should be linked to the boundaries conditions class
    std::vector<double> C_n(size_matrix); //Value taken from boundaries conditions
    std::vector<double> C_n_next(size_matrix); //Value taken from boundaries conditions
    
    //construct the vector X to be determined by the system of equations, by setting to the value Final payoff of the option --> should be linked to the payoff class
    std::vector<double> X(size_matrix);
    // X = Value taken from final payoff class
    
    //Find the value of X_(n) by solving the system: AX_(n) = BX_(n+1) + C_(n+1) - C_(n)
    X = xt::linalg::solve(A, xt::linalg::dot(B, X) + C_n_next - C_n );
        
        
};
        
        
//all in function to get the price of the option and the greeks ?
std::vector<double> resolution(){};

//function to compute the greeks ?
