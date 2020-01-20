#include "Resolution.hpp"

Solve::Resolution(double dt, double dx, long nb_step_time, double v, double r, double theta, long nb_step_spot, std::vector<double> conditions, std::vector<double> m_init_vector){

    
//Converts the matrix of boundaries conditions into an xt::array type
xt::xarray<double> C_n = xt::adapt(conditions);
        
//Converts the vector X_temp into an xt::array type vector named X
xt::xarray<double> X = xt::adapt(m_init_vector);
    

for (int t = 2; t <= nb_step_time; ++t){
    
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
        B_X.reshape({nb_step_spot,1});
    
        //Determines the value of X_(n) by solving the system: AX_(n) = B.X_(n+1) + C_(n+1) - C_(n)
        X = xt::linalg::solve(A, B_X + xt::view(C_n, xt::all(), xt::range(nb_step_time-t+1,nb_step_time-t+2)) - xt::view(C_n, xt::all(), xt::range(nb_step_time-t,nb_step_time-t+1)));
};
}
