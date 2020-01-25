
#include "Resolve.hpp"
#include <vector>
#include <cmath>
#include <limits>


 solver::solver(mesh grid, std::vector<std::vector<double>> res):m_mesh(grid),m_results(res){};
 
 void solver::solve_X(mesh grid, double theta, const std::vector<double>& boundaries, std::vector<std::vector<double>> vol_mat,std::vector<std::vector<double>> rate_mat)
 {
	 
	 
	 std::vector<std::vector<double>> results;
	 std::vector<std::vector<double>> vol(vol_mat);
	 std::vector<std::vector<double>> rate(rate_mat);
	 long T = grid.Getvector_time().size(); 
	 std::size_t m = grid.Getvector_stock().size();
	 std::vector<double> init_f = grid.get_init_vector(); // get the initial X_T
	 
	 
	 
	 init_f.erase(init_f.begin());
	 init_f.pop_back();
	 
	 std::vector<std::vector<double>> matrix_boundaries = transform_matrix(boundaries,m-2); //get the boundaries matrix 
	 std::vector<double> cond_diff(m-2);
	 std::vector<double> solution_vector(m-2, 0.0); //init of the solution vector 
	 std::vector<double> up_A;
	 std::vector<double> low_A;
	 std::vector<double> mid_A;
	 std::vector<double> up_B;
	 std::vector<double> low_B;
	 std::vector<double> mid_B;
	 
	 std::cout << "price at time T" << std::endl;
	 
	 print(init_f);
	 
	 solver solver(grid, results);
	 
	 double BetaN = 1;
	 double BetaN1 = 1;
	 double AlphaN = 1;
	 double AlphaN1 = 1;
	 
	 	 //initialisation de la matrice B à maturité 
	 up_B = solver.Upper_diag_coeff(grid,false,theta,vol.back(), rate.back());
	 low_B = solver.Lower_diag_coeff(grid,false,theta,vol.back(), rate.back());
	 mid_B = solver.Mid_diag_coeff(grid,false,theta,vol.back(), rate.back());
	 
	up_A = solver.Upper_diag_coeff(grid,true,theta,vol.back(), rate.back());
	low_A = solver.Lower_diag_coeff(grid,true,theta,vol.back(), rate.back());
	mid_A = solver.Mid_diag_coeff(grid,true,theta,vol.back(), rate.back());
	
	std::cout << "test" << std::endl;
	 
	 for(size_t s = 0; s<m-2; ++s)
	 {
		 
		 if(s==0){
			 
			 BetaN = low_A[0];

			 BetaN1 = low_B[0];
			 
			 cond_diff[s] = BetaN1*matrix_boundaries[s].back()-BetaN*matrix_boundaries[s][T-1];
		 }
		 else if( s ==m-3){
			 
			 AlphaN = up_A.back();
			 AlphaN1 = up_B.back();
			 
			 cond_diff[s] = AlphaN1*matrix_boundaries[s].back()-AlphaN*matrix_boundaries[s][T-1];
			 
		 }
		 
		 //cond_diff[s] = matrix_boundaries[s].back()-matrix_boundaries[s][T-1];
		 
		 else{cond_diff[s] =0;};
	 }; //initialisation du vecteur de C(n+1) - C(n) 
	 
	 

	 

	std::vector<double> B = solver.BX_vector(up_B,mid_B,low_B,cond_diff,init_f);

	rate.pop_back();
	vol.pop_back();
	 
	 
	 for(int i = grid.Getvector_time().size() - 1; i != 0; --i)
	 {
		 
		
		up_A = solver.Upper_diag_coeff(grid,true,theta,vol.back(), rate.back());
		low_A = solver.Lower_diag_coeff(grid,true,theta,vol.back(), rate.back());
		mid_A = solver.Mid_diag_coeff(grid,true,theta,vol.back(), rate.back());
			
		solver.thomas_algorithm(up_A, mid_A, low_A, B, solution_vector);
		
		results.push_back(solution_vector); //create a matrix containing all the prices computed by the solver in order to be displayed afterwards
		
		
		for(size_t s = 0; s<m-2; ++s){
		 	
		if(s==0){
			 
			 BetaN = low_A[0];

			 BetaN1 = low_B[0];
			 
			 cond_diff[s] = BetaN1*matrix_boundaries[s][i]-BetaN*matrix_boundaries[s][i-1];
		 }
		 else if( s==m-3){
			 
			 AlphaN = up_A.back();
			 AlphaN1 = up_B.back();
			 
			 cond_diff[s] = AlphaN1*matrix_boundaries[s][i]-AlphaN*matrix_boundaries[s][i-1];
			 
		 }
		 //cond_diff[s] = matrix_boundaries[s][i]-matrix_boundaries[s][i-1];
		 else{ cond_diff[s] =0;};
			 }; 
			 
			 
		//init_f = solution_vector;
		
		std::cout << "price at time " << i-1 << " is " << std::endl;
		
		print(solution_vector);
		
		up_B = solver.Upper_diag_coeff(grid, false,theta,vol.back(), rate.back());
	    low_B = solver.Lower_diag_coeff(grid, false,theta,vol.back(), rate.back());
		mid_B = solver.Mid_diag_coeff(grid, false,theta,vol.back(), rate.back());
	 
	 
	    B = solver.BX_vector(up_B,mid_B,low_B,cond_diff,solution_vector);
		
		rate.pop_back();
		vol.pop_back();
		
		 
	 }; 
	  
 }; 
 
std::vector<std::vector<double>> solver::get_price(){
	
	return m_results;
};
 
 std::vector<double> solver::BX_vector(std::vector<double> upper, std::vector<double> mid, std::vector<double> low, std::vector<double> bound_diff, std::vector<double> Fn1){
//this procedure creates the right_hand vector of the AX = D equation based on BX(n+1) + C(n+1) - C(n)
	std::size_t N = mid.size(); // as the resolution si between f1 and fN-1
	
	std::vector<double> BX;
	
	
	BX.push_back(mid[0]*Fn1[0] + upper[0]*Fn1[1] +bound_diff[0]);
	
	
	for(size_t i = 1; i < N-1; i++){

   	BX.push_back(mid[i]*Fn1[i] + upper[i]*Fn1[i+1] +low[i-1]*Fn1[i-1] + bound_diff[i]);
	
	};
	
	BX.push_back(mid.back()*Fn1.back() + low.back()*Fn1[N-1] + bound_diff.back());
	
	return BX;
};
	
//set of function to define the A matrix (3 better than just one huge matrix ?) 
std::vector<double> solver::Mid_diag_coeff(mesh grid, bool A,double theta, std::vector<double> sigma, std::vector<double> rate){
	
	sigma.pop_back();
	sigma.erase(sigma.begin());
	
	rate.pop_back();
	rate.erase(rate.begin());
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	//std::vector<double> sigma_init = param.get_Vol(); //to get the volatility 
	//std::vector<double> rate_init = param.get_Rate(); //the get the rate 
	//double theta = param.get_Theta(); //to get the theta 
	double size_gamma = grid.Getvector_stock().size()-2; //no minus-1 as we are on diagonal 
	
	
	
		//create the vector that holds the diagonal 
	std::vector<double> gamma_coefficient(size_gamma);
	
	for (std::size_t i = 0; i < size_gamma; ++i){
		
		gamma_coefficient[i] =  std::pow(sigma[i],2)/std::pow(dx,2) + rate[i];
		
		if (A==false){
			
			gamma_coefficient[i] = -dt*(1-theta)*gamma_coefficient[i] + 1;
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		}
		
		else{
			
			gamma_coefficient[i] = dt*theta*gamma_coefficient[i] + 1;
		};
	};
	
	return gamma_coefficient;
	
}; 


std::vector<double> solver::Upper_diag_coeff(mesh grid, bool A,double theta,std::vector<double> sigma, std::vector<double> rate){
	
	
	sigma.pop_back();
	sigma.pop_back();
	sigma.erase(sigma.begin());
	
	rate.pop_back();
	rate.pop_back();
	rate.erase(rate.begin());
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	//std::vector<double> sigma_init = param.get_Vol(); //to get the volatility 
	//std::vector<double> rate_init = param.get_Rate(); //the get the rate 
	//double theta = param.get_Theta(); //to get the theta 
	double size_alpha = grid.Getvector_stock().size()-3; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> alpha_coefficient(size_alpha);
	
	for (std::size_t i = 0; i < size_alpha; ++i){
		
		alpha_coefficient[i] =  (-std::pow(sigma[i],2))/(2*std::pow(dx,2)) + std::pow(sigma[i],2)/(4*std::pow(dx,2)) - (rate[i])/(2*dx);
		
		 
		if (A==false){
			
			alpha_coefficient[i] = -dt*(1-theta)*alpha_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		}
		
		else{
			
			alpha_coefficient[i] = dt*theta*alpha_coefficient[i];
		};
	};
	
	return alpha_coefficient;
	
};


std::vector<double> solver::Lower_diag_coeff(mesh grid, bool A,double theta,std::vector<double> sigma, std::vector<double> rate){
	
	
	sigma.pop_back();
	sigma.erase(sigma.begin());
	sigma.erase(sigma.begin());
	
	rate.pop_back();
	rate.erase(rate.begin());
	rate.erase(rate.begin());
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	//std::vector<double> sigma_init = param.get_Vol(); //to get the volatility 
	//std::vector<double> rate_init = param.get_Rate(); //the get the rate 
	//double theta = param.get_Theta(); //to get the theta 
	double size_beta = grid.Getvector_stock().size()-3; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> beta_coefficient(size_beta);
	
	for (std::size_t i = 0; i < size_beta; ++i){
		
		beta_coefficient[i] =  (-std::pow(sigma[i],2))/(2*std::pow(dx,2)) + std::pow(sigma[i],2)/(4*std::pow(dx,2)) + (rate[i])/(2*dx);
		
		
		if (A==false){
			
			beta_coefficient[i] = -dt*(1-theta)*beta_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		}
		
		else{
			
			beta_coefficient[i] = dt*theta*beta_coefficient[i];
		};
	};
	
	return beta_coefficient;
	
};

//this is the thomas algo for inverting the matrix  

void solver::thomas_algorithm(const std::vector<double>& upper_diag, const std::vector<double>& mid_diag, const std::vector<double>& lower_diag, const std::vector<double>& f_n1, std::vector<double>& f_sol) {
  size_t nb_spot = f_n1.size();
  
  std::vector<double>  lower_diag2(lower_diag);
  std::vector<double>  upper_diag2(upper_diag);
  
  upper_diag2.push_back(0.0);
  lower_diag2.insert(lower_diag2.begin(),0.0);
                                                                                                                                                                                                                                                                                                                                                 
  std::vector<double> c_star(nb_spot, 0.0);
  std::vector<double> d_star(nb_spot, 0.0);                                                                                                                                                    
  c_star[0] = lower_diag2[0] / mid_diag[0];
  d_star[0] = f_n1[0] / mid_diag[0];

  //forward sweep                                                                                                                                                  
  for (int i=1; i<nb_spot; i++) {
    double m = 1.0 / (mid_diag[i] - upper_diag2[i] * c_star[i-1]);
    c_star[i] = lower_diag2[i] * m;
    d_star[i] = (f_n1[i] - upper_diag2[i] * d_star[i-1]) * m;
  }
  
  f_sol.back() = (d_star.back() - lower_diag2.back()*d_star[nb_spot-2]) / (mid_diag.back() - lower_diag2.back()*upper_diag2.back());

  //reverse sweep, used to update the solution vector f                                                                                                                                                 
  for (int i=nb_spot-1; i-- > 0; ) {
    f_sol[i] = d_star[i] - c_star[i] * f_n1[i+1];
  }
}; 


solver::~solver(){}; 