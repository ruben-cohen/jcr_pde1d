#include <fstream>
#include "fdm.h"

//implement a constructor for the abstract base class FDMBase. This is because it is actually storing member data and so we need to call this constructor from a derived class' member initialisation list
FDMBase::FDMBase(double _x_dom, unsigned long _J,
                 double _t_dom, unsigned long _N,
                 PDE* _pde) 
  : x_dom(_x_dom), J(_J), t_dom(_t_dom), N(_N), pde(_pde) {}
  
  
//The constructor for FDMEulerExplicit  calls the methods to fill in the step sizes and initial conditions
FDMEulerExplicit::FDMEulerExplicit(double _x_dom, unsigned long _J,
                                   double _t_dom, unsigned long _N,
                                   PDE* _pde) 
  : FDMBase(_x_dom, _J, _t_dom, _N, _pde) 
{
  calculate_step_sizes();
  set_initial_conditions();
}
//For N temporal discretisation points we have N−1 intervals. Similarly for the spatial discretisation. This method calculates these steps values for later use
void FDMEulerExplicit::calculate_step_sizes() 
{
  dx = x_dom/static_cast<double>(J-1);
  dt = t_dom/static_cast<double>(N-1);
}
//Now that the step sizes are set it is time to pre-fill the initial condition.
//All of the spatial arrays are set to have J points and are zeroed. Then the method loops these arrays and uses the pointer to the PDE to obtain the initial condition as a function of spot. 
//We also have a useful "helper" array, x_values which stores the spot value at each discretisation point to save us calculating it every time step.
//Finally, we set the current and previous times to zero:
void FDMEulerExplicit::set_initial_conditions() 
{
  // Spatial settings
  double cur_spot = 0.0;

  old_result.resize(J, 0.0);
  new_result.resize(J, 0.0);
  x_values.resize(J, 0.0);

  for (unsigned long j=0; j<J; j++) 
  {
    cur_spot = static_cast<double>(j)*dx;
    old_result[j] = pde->init_cond(cur_spot);
    x_values[j] = cur_spot;
  }

  // Temporal settings
  prev_t = 0.0;
  cur_t = 0.0;
}
//Now that the initial conditions are set we can begin the time-marching.
// However, in each time step we must first set the boundary conditions. 
//In this instance the edge points (at index 0 and index J−1) are set using the PDE boundary_left and boundary_right method. This is an example of a Dirichlet condition.
// More sophisticated boundary conditions approximate the derivative at these points (Neumann conditions), although we won't be utilising these conditions here:
void FDMEulerExplicit::calculate_boundary_conditions() 
{
  new_result[0] = pde->boundary_left(prev_t, x_values[0]);
  new_result[J-1] = pde->boundary_right(prev_t, x_values[J-1]);
}
//Main method to solve the inner domain of our mesh
// A loop is carried out over the spatial cells - but only those away from the boundary,
//The next step is to calculate the α, β and γ coefficients which represent the algebraic rearrangement of the derivative discretisation.
//Notice that we need to obtain certain coefficients via pointer dereferencing of the underlying PDE
//Once α, β and γ are defined we can use the finite differencing to update the solution into the new_result vector. This formula will also be clear from the mathematical algorithm above:
void FDMEulerExplicit::calculate_inner_domain() 
{
  // Only use inner result indices (1 to J-2)
  for (unsigned long j=1; j<J-1; j++) {
    // Temporary variables used throughout
    double dt_sig = dt * (pde->diff_coeff(prev_t, x_values[j]));
    double dt_sig_2 = dt * dx * 0.5 * (pde->conv_coeff(prev_t, x_values[j]));

    // Differencing coefficients (see \alpha, \beta and \gamma in text)
    alpha = dt_sig - dt_sig_2;
    beta = dx * dx - (2.0 * dt_sig) + (dt * dx * dx * (pde->zero_coeff(prev_t, x_values[j])));
    gamma = dt_sig + dt_sig_2;

    // Update inner values of spatial discretisation grid (Explicit Euler)
    new_result[j] = ( (alpha * old_result[j-1]) + 
                      (beta * old_result[j]) + 
                      (gamma * old_result[j+1]) )/(dx*dx) - 
      (dt*(pde->source_coeff(prev_t, x_values[j])));
  }
}

//Now that all of the prior pure virtual methods have been given an implementation in the derived FDMEulerExplicit class it is possible to write the step_march method to "wrap everything together".

void FDMEulerExplicit::step_march() 
{ 
  std::ofstream fdm_out("fdm.csv");

  while(cur_t < t_dom) {
    cur_t = prev_t + dt;
    calculate_boundary_conditions();
    calculate_inner_domain();
    for (int j=0; j<J; j++) {
      fdm_out << x_values[j] << " " << prev_t << " " << new_result[j] << std::endl;
    }
    
    old_result = new_result;
    prev_t = cur_t;
  }

  fdm_out.close();
}
