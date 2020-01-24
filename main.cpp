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
	double r = 0.0;   // Risk-free rate (5%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 1.00;    // One year until expiry
	double theta_ = 0.5;
	
	// mesh discretisation parameters
	long nb_step_spot =50;    // Spot goes from [0.0, 1.0]
	long nb_step_time = 10; 
	
	
	// Create the PayOff object (strike as parameter)
	PayOff* option = new PayOffCall(K);
	// Create the mesh object (parameters (S,Vol, Nb time step, nb stock step)
	mesh grille(S,T,v,nb_step_time,nb_step_spot,option);
	
	//Create PDE object (to solve the problem, we provide to the class, the grid previoulsy defined in the class PDE) 
	//PDE m(option,grille.getdx(),grille.getdt(),grille.Getvector_time(),grille.Getvector_stock());
	
	//test of the boundaries : 
	//Parameters par(v, r, theta_);
	long s = grille.Getvector_stock().size();
	long _t_ = grille.Getvector_time().size();
	std::vector<double> sigma(s,v);
	std::vector<double> rate(s,r);
	std::vector<double> K_v{0.1,0.1,0.1,0.1};
	Derichtlet c(grille, rate);
	Neumann c2(grille, theta_, sigma, rate,K_v);
	
	//std::cout<< "size is" << c2.get_cond().size() << std::endl;
	//print(c2.get_cond());
	
	
	std::vector<std::vector<double>> test = transform_matrix(c.get_cond(),s);
	
	
	
/* 		for(long i=0; i<test.size()-1;i++){
		
			for(long j=0; j<test[j].size()-1;j++){
		
					std::cout << test[j][i]  << std::endl;
					//std::cout << rate_mat.back()[i] << std::endl;
		
				};
				
		
		
	}; */
	
	
	
	std::vector<std::vector<double>> vol_mat;
	std::vector<std::vector<double>> rate_mat;
	
	for(long i=0; i<_t_;i++){
		
		vol_mat.push_back(sigma);
		rate_mat.push_back(rate);
	};
	
/* 	for(long i=0; i<vol_mat.size()-1;i++){
		
			for(long j=0; j<vol_mat[j].size()-1;j++){
		
					std::cout << vol_mat.back()[i] << std::endl;
					std::cout << rate_mat.back()[i] << std::endl;
		
				};
		
	}; */
	std::vector<double> up_B;
	std::vector<double> mid_B;
	std::vector<double> low_B;
	std::vector<double> B;
	
	std::vector<double> init_f(grille.get_init_vector());
	
	std::vector<std::vector<double>> res;
	
	solver sol(grille, res);
	

	up_B = sol.Upper_diag_coeff(grille,false,theta_,sigma, rate);
    low_B = sol.Lower_diag_coeff(grille,false,theta_,sigma, rate);
	mid_B = sol.Mid_diag_coeff(grille,false,theta_,sigma, rate);

	 
	std::vector<double> cond_test;
	
	for(int i =0; i <test.size();i++){
		
		
		cond_test.push_back(test[i].back());
	}
	
	//std::cout << "size_vector test " << cond_test.size() <<std::endl;
	
	//print(init_f);
	
	init_f.erase(init_f.begin());
	init_f.pop_back();
	 
	B = sol.BX_vector(up_B,mid_B,low_B,cond_test,init_f);
	
	//print(B);
	
	//std::cout << B.size() << std::endl;
	
/* 	std::cout << "B" << std::endl;
	

	//print(mid_B);
	std::cout << mid_B.size() << std::endl;
	//print(up_B);
	
	std::cout << up_B.size() << std::endl;
	//print(low_B);

	std::cout << low_B.size() << std::endl; */
	
	
	
    sol.solve_X(grille, theta_,c.get_cond(), vol_mat, rate_mat);
	
	
	//print(K_v);
	
	//std::vector<double> m_data(2 * 3);
	
	//m_data
	// m_data(1
	
	
	//xt::xarray<double> diago = xt::eye(20,-1);
	// xt::xtensor<double, 1> K_v = {0.1, 0.1, 0.1, 0.1 };
	// xt::xtensor<double, 1> K_v2 = {0.1, 0.1, 0.1, 0.1 };
	// xt::xarray<double> res = xt::linalg::dot(K_v, K_v2);
	//std::vector<double> der = c.get_cond();
	// std::cout<< res()<< std::endl;
	//xt::xarray<double> res = xt::linalg::dot(K_v, K_v2);
	
	//Solve sv(grille, par, der);
	
	//xt::xarray<double> price = sv.Get_FX_n();
	
	//price.reshape({grille.Getvector_stock().size(),1});
	//std::cout<< diago << std::endl;
	//std::cout<< "solver" << std::endl;
	//std::cout << price << std::endl;
	
/* 	std::cout<< "fonction calcul cond init:" << std::endl;
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
	print(c2.get_cond());
	std::cout<< "Size vecteur dirichlet:" << std::endl;
	std::cout<< c2.get_cond().size() <<std::endl;
	
	
	std::cout<< "vecteur Neumann:" << std::endl;
	print(c2.get_cond()); */

	delete option;



	return 0;
	
}
