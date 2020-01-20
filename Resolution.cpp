#include "Resolution.hpp"

Solve::Solve(mesh _grid, Parameters _param, bound_conditions _conditions){

double dt = _grid.getdt();
double dx = _grid.getdx();
std::vector<double> m_init_vector = _grid.get_init_vector();
double v = _param.Get_Vol();
double r = _param.Get_Rate();
double theta = _param.Get_Theta();
long nb_step_time = _grid.Getvector_time().size();
long nb_step_stock = _grod.Getvector_stock().size();


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
