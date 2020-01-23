#include <iostream>
#include "closed_form.hpp"
#include <vector>
//#include "Resolution.hpp"
#include <cmath>
//Allows to create and manipulate xarrays
// #include "xtensor/xarray.hpp"
// #include "xtensor/xio.hpp"
// #include "xtensor/xview.hpp"
// #include "xtensor/xadapt.hpp"
// #include "xtensor/xeval.hpp"
// // Allows to perform linear algebra operations on xarrays
// #include "xtensor-blas/xlinalg.hpp"


////
int main(int argc, char* argv[])
{
	// Create the option parameters
	double S = 100;
	double K = 100;  // Strike price
	double r = 0.05;   // Risk-free rate (5%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 1.00;    // One year until expiry
	double theta_ = 0.5;
	
	// mesh discretisation parameters
	long nb_step_spot =20;    // Spot goes from [0.0, 1.0]
	long nb_step_time = 20; 
	
	
	// Create the PayOff object (strike as parameter)
	PayOff* option = new PayOffCall(K);
	// Create the mesh object (parameters (S,Vol, Nb time step, nb stock step)
	mesh grille(S,T,v,nb_step_time,nb_step_spot,option);
	
	//Create PDE object (to solve the problem, we provide to the class, the grid previoulsy defined in the class PDE) 
	//PDE m(option,grille.getdx(),grille.getdt(),grille.Getvector_time(),grille.Getvector_stock());
	
	//test of the boundaries : 
	//Parameters par(v, r, theta_);
	//std::vector<double> K_v{0.1,0.1,0.1,0.1};
	//std::vector<double> K_v2{0.1,0.1,0.1,0.1};
	//Derichtlet c(grille,par);
	//Neumann c2(grille,par,K_v);
	//print(K_v);
	
	//std::vector<double> m_data(2 * 3);
	
	//m_data
	// m_data(1
	
	
	//xt::xarray<double> diago = xt::eye(20,-1);
	// xt::xtensor<double, 1> K_v = {0.1, 0.1, 0.1, 0.1 };
	// xt::xtensor<double, 1> K_v2 = {0.1, 0.1, 0.1, 0.1 };
	// xt::xarray<double> res = xt::linalg::dot(K_v, K_v2);
	std::vector<double> der = c.get_cond();
	// std::cout<< res()<< std::endl;
	//xt::xarray<double> res = xt::linalg::dot(K_v, K_v2);
	
	//Solve sv(grille, par, der);
	
	//xt::xarray<double> price = sv.Get_FX_n();
	
	//price.reshape({grille.Getvector_stock().size(),1});
	//std::cout<< diago << std::endl;
	//std::cout<< "solver" << std::endl;
	//std::cout << price << std::endl;
	
	std::cout<< "fonction calcul cond init:" << std::endl;
	std::cout << grille.init_cond(S) << std::endl;

	std::cout<< "Vecteur de prix (log):" << std::endl;
	print(grille.Getvector_stock());
	std::cout<< "Size vecteur stock:" << std::endl;
	std::cout<< grille.Getvector_stock().size() <<std::endl;
	
	std::cout<< "Vecteur de temps:" << std::endl;
	print(grille.Getvector_time());
	std::cout<< "Size vecteur temps:" << std::endl;
	std::cout<< grille.Getvector_time().size() <<std::endl;
	
	std::cout<< "Vecteur cond init (not log):" << std::endl;
	print(grille.get_init_vector());
	std::cout<< "Size vecteur cond ini:" << std::endl;
	std::cout<< grille.get_init_vector().size() <<std::endl;
	
	std::cout<< "dx:" << std::endl;
	std::cout << grille.getdx() << std::endl;
	
	std::cout<< "dt:" << std::endl;
	std::cout << grille.getdt() << std::endl;
	
	std::cout<< "vecteur dirichlet:" << std::endl;
	print(c.get_cond());
	std::cout<< "Size vecteur dirichlet:" << std::endl;
	std::cout<< c.get_cond().size() <<std::endl;
	
	
	// std::cout<< "vecteur Neumann:" << std::endl;
	// print(c2.get_cond());

	delete option;



	return 0;
	
}
