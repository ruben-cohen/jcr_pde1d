
#include <vector>

// Finite Difference Method - Abstract Base Class
class FDMBase 
{
	 protected:
	  PDE* pde;

	  // Space discretisation
	  double x_dom_max;      // Spatial extent [0.0, x_dom]
	  double x_dom_min;      // Spatial extent borne min
	  unsigned long J;   // variable pour le nombre de points
	  double dx;         // Step size spatial
	  std::vector<double> x_values;  // vecteur pour les coordonnées dimension spatiale

	  // Time discretisation
	  double t_dom;      // dimension temps [0.0, t_dom]
	  unsigned long N;   //  nombre de points sur la dimension temps
	  double dt;         // Temporal step size (calculated from above)

	  // Time-marching
	  double prev_t, cur_t;   // temps présent et temps passé

	  // coefficients de notre matrice à inverser
	  double alpha, beta, gamma;

	  // Storage
	  std::vector<double> new_result;   // New solution (becomes N+1)
	  std::vector<double> old_result;   // Old solution (becomes N)

	  // Constructor
	  FDMBase(double _x_dom, unsigned long _J,
			  double _t_dom, unsigned long _N,
			  PDE* _pde);

	  // Override these virtual methods in derived classes for 
	  // specific FDM techniques, such as explicit Euler, Crank-Nicolson, etc.
	  virtual void calculate_step_sizes() = 0;
	  virtual void set_initial_conditions() = 0;
	  virtual void calculate_boundary_conditions() = 0;
	  virtual void calculate_inner_domain() = 0;

	 public:
	  // Carry out the actual time-stepping
	  virtual void step_march() = 0;
};

class FDMEulerExplicit : public FDMBase 
	
{
	 protected:
	  void calculate_step_sizes();
	  void set_initial_conditions();
	  void calculate_boundary_conditions();
	  void calculate_inner_domain();

	 public:
	  FDMEulerExplicit(double _x_dom, unsigned long _J,
					   double _t_dom, unsigned long _N,
					   ConvectionDiffusionPDE* _pde);

	  void step_march();
};
