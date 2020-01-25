#include <cmath>
#include <limits>
#include <algorithm>
#include "Greeks.hpp"


Greeks::Greeks(mesh grid, solver, Neumann, Derichtlet){

double dt = grid.getdt(); //get the time step
double dx = grid.getdx(); //get the stock step
long T = grid.Getvector_time().size(); //get the number of time steps
double nb_spot_values = grid.Getvector_stock().size(); //get the number of spot values
std::vector<double> price_matrix = solver.get_price(); //get the matrix of option prices
std::vector<double> Derichtlet.get_cond() //get the Dirichlet condition vector
std::vector<double> Neumann.get_cond() //get the Neumann condition vector
    
std::vector<double> delta;
std::vector<double> gamma;
std::vector<double> theta;

// Boundaries conditions
// Assuming that if the user imputs a Dirichlet condition, then the Neumann conditions (K1...) are null
delta[0] = get_K3()[0];
gamma[0] = get_K4()[0];
theta[0] = (price_matrix[1][0] - price_matrix[0][0])/dt;

delta[nb_spot_values] = get_K1()[0];
gamma[nb_spot_values] = get_K2()[0];
theta[nb_spot_values] = (price_matrix[1][nb_spot_values] - price_matrix[0][nb_spot_values])/dt;
    
    
 for(size_t i = 1; i < nb_spot_values-1; ++i){
 
     delta[i] = (price_matrix[0][i+1] - price_matrix[0][i-1])/(2*dx);
     gamma[i] = (price_matrix[0][i+1] - 2*price_matrix[0][i] + price_matrix[0][i+1])/(dx*dx);
     theta[i] = (price_matrix[1][i] - price_matrix[0][i])/dt;
     
 };

};

std::vector<double> Greeks::get_delta(){
    
    return delta;
};

std::vector<double> Greeks::get_gamma(){
    
    return gamma;
};

std::vector<double> Greeks::get_theta(){
    
    return theta;
};
