
#include 'Resolve.hpp'
#include <vector>
#include <cmath>


 solver::solver(mesh grid, double theta, const std::vector<double>& boundaries, std::vector<std::vector<double>> vol_mat,std::vector<std::vector<double>> rate_mat)
	m_mesh(grid)
 {
	 std::vector<std::vector<double>> results;
	 std::vector<std::vector<double>> vol(vol_mat);
	 std::vector<std::vector<double>> rate(rate_mat);
	 long T = grid.Getvector_time().size() 
	 std::size_t m = grid.Getvector_stock().size();
	 std::vector<double> init_f = grid.get_init_vector(); // get the initial X_T
	 std::vector<std::vector<double>> matrix_boundaries = transform_matrix(boundaries); //get the boundaries matrix 
	 std::vector<double> cond_diff(m);
	 std::vector<double> solution_vector(m-2, 0.0); //init of the solution vector 
	 std::vector<double> up_A;
	 std::vector<double> low_A;
	 std::vector<double> mid_A;
	 std::vector<double> up_B;
	 std::vector<double> low_B;
	 std::vector<double> mid_B;
	 
	 
	 for(size_t s = 0; s<m; ++s){
		 
		 cond_diff[s] = matrix_boundaries[s][T]-matrix_boundaries[s][T-1];
	 }; //initialisation du vecteur de C(n+1) - C(n) 
	 
	 //initialisation de la matrice B à maturité 
	 up_B = Upper_diag_coeff(grid, param,false,param.Gettheta(),vol.back(), rate.back());
	 low_B = Lower_diag_coeff(grid, param,false,param.Gettheta(),vol.back(), rate.back());
	 mid_B = Mid_diag_coeff(grid, param,false,param.Gettheta(),vol.back(), rate.back());
	 
	 
	std::vector<double> B = BX(up,mid,low,cond_diff,init_f);

	rate.pop_back();
	vol.pop_back();
	 
	 
	 for(int i = grid.Getvector_time().size() - 1; i > 0; --i)
	 {
		 
		
		up_A = Upper_diag_coeff(grid, param,true,theta,vol.back(), rate.back());
		low_A = Lower_diag_coeff(grid, param,true,theta,vol.back(), rate.back());
		mid_A = Mid_diag_coeff(grid, param,true,theta,vol.back(), rate.back());
			
		thomas_algorithm(up_A, mid_A, low_A, B, solution_vector);
		results.push_back(solution_vector); //create a matrix containing all the prices computed by the solver in order to be displayed afterwards
		
		for(size_t s = 0; s<m; ++s){
		 
				cond_diff[s] = matrix_boundaries[s][i]-matrix_boundaries[s][i-1];
			 }; 
			 
		init_f = solution_vector;
		
		up_B = Upper_diag_coeff(grid, param,false,theta,vol.back(), rate.back());
	    low_B = Lower_diag_coeff(grid, param,false,theta,vol.back(), rate.back());
		mid_B = Mid_diag_coeff(grid, param,false,theta,vol.back(), rate.back());
	 
	 
	    B = BX(up,mid,low,cond_diff,init_f);
		
		rate.pop_back();
		vol.pop_back();
		
		 
	 }; 
	 
	 
	 
	 
 }; 
 
const std::vector<std::vector<double>> get_price(){
	
	return results;
};
 
 std::vector<double> BX_vector(std::vector<double> upper, std::vector<double> mid, std::vector<double> low, std::vector<double> bound_diff, std::vector<double> Fn1){
//this procedure creates the right_hand vector of the AX = D equation based on BX(n+1) + C(n+1) - C(n)
	std::size_t N = mid.size(); // as the resolution si between f1 and fN-1
	
	std::vector<double> BX(N);
	
	for(size_t i = 0; i < N; i++){
		
		
		switch (i) {
        case 0: BX[i] = mid[0]*Fn1[0] + upper[0]*Fn1[0] +bound_diff[0];
		case N: BX[i] =  mid.back()*Fn1.back() + lower.back()*Fn1[N-1] + bound_diff.back();
        default: BX[i] = mid[i]*Fn1[i] + upper[i]*Fn1[i+1] +lower[i-1]*Fn1[i-1] + bound_diff(i);
    };
	
	return BX;
		
	};
};
	
//set of function to define the A matrix (3 better than just one huge matrix ?) 
std::vector<double> solver::Mid_diag_coeff(mesh grid, double theta,bool A,double theta, std::vector<double> sigma, std::vector<double> rate){
	
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	//std::vector<double> sigma_init = param.get_Vol(); //to get the volatility 
	//std::vector<double> rate_init = param.get_Rate(); //the get the rate 
	//double theta = param.get_Theta(); //to get the theta 
	double size_gamma = Getvector_stock().size(); //no minus-1 as we are on diagonal 
	
	
	
		//create the vector that holds the diagonal 
	std::vector<double> gamma_coefficient(size_gamma);
	
	for (std::size_t i = 0; i < size_gamma; ++i){
		
		gamma_coefficient[i] =  std::pow(sigma[i],2)/std::pow(dx,2) + rate[i];
		
		if (A==false){
			
			gamma_coefficient[i] = -dt*(1-theta)*gamma_coefficient[i] + 1;
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		};
		
		else{
			
			gamma_coefficient[i] = dt*theta*alpha_coefficient[i] + 1;
		};
	};
	
	return gamma_coefficient;
	
}; 


std::vector<double> solver::Upper_diag_coeff(mesh grid, double theta,bool A,double theta,std::vector<double> sigma, std::vector<double> rate){
	
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	//std::vector<double> sigma_init = param.get_Vol(); //to get the volatility 
	//std::vector<double> rate_init = param.get_Rate(); //the get the rate 
	//double theta = param.get_Theta(); //to get the theta 
	double size_alpha = Getvector_stock().size()-1; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> alpha_coefficient(size_alpha);
	
	for (std::size_t i = 0; i < size_alpha; ++i){
		
		alpha_coefficient[i] =  (-std::pow(sigma[i],2))/(2*std::pow(dx,2)) + std::pow(sigma[i],2)/(4*std::pow(dx,2)) + (-rate[i])/(2*dx);
		
		 
		if (A==false){
			
			alpha_coefficient[i] = -dt*(1-theta)*alpha_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		};
		
		else{
			
			alpha_coefficient[i] = dt*theta*alpha_coefficient[i];
		};
	};
	
	return alpha_coefficient;
	
};


std::vector<double> solver::Lower_diag_coeff(mesh grid, double theta,bool A,double theta,std::vector<double> sigma, std::vector<double> rate){
	
	
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	//std::vector<double> sigma_init = param.get_Vol(); //to get the volatility 
	//std::vector<double> rate_init = param.get_Rate(); //the get the rate 
	//double theta = param.get_Theta(); //to get the theta 
	double size_beta = Getvector_stock().size()-1; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> beta_coefficient(size_beta);
	
	for (std::size_t i = 0; i < size_beta; ++i){
		
		beta_coefficient[i] =  (-std::pow(sigma[i],2))/(2*std::pow(dx,2)) + std::pow(sigma[i],2)/(4*std::pow(dx,2)) + (rate[i])/(2*dx);
		
		
		if (A==false){
			
			beta_coefficient[i] = -dt*(1-theta)*beta_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		};
		
		else{
			
			beta_coefficient[i] = dt*theta*beta_coefficient[i];
		};
	};
	
	return beta_coefficient;
	
};

//this is the thomas algo for inverting the matrix  

void solver::thomas_algorithm(const std::vector<double>& upper_diag,
                      const std::vector<double>& mid_diag,
                      const std::vector<double>& lower_diag,
                      const std::vector<double>& f_n1,
                      std::vector<double>& f_sol) {
  size_t nb_spot = f_n1.size();

  // Create the temprary vectors to store new coef                                                                                                                                                                                                                                                                                                                                                
  std::vector<double> c_star(nb_spot, 0.0);
  std::vector<double> d_star(nb_spot, 0.0);                                                                                                                                                    
  c_star[0] = lower_diag[0] / mid_diag[0];
  d_star[0] = f_n1[0] / mid_diag[0];

  //forward sweep                                                                                                                                                  
  for (int i=1; i<nb_spot; i++) {
    double m = 1.0 / (mid_diag[i] - upper_diag[i] * c_star[i-1]);
    c_star[i] = lower_diag[i] * m;
    d_star[i] = (f_n1[i] - upper_diag[i] * d_star[i-1]) * m;
  }

  //reverse sweep, used to update the solution vector f                                                                                                                                                 
  for (int i=nb_spot-1; i-- > 0; ) {
    f_sol[i] = d_star[i] - c_star[i] * f_n1[i+1];
  }
}; //this is a void function because it modifies the price 

std::vector<double> resolution(){};
	
//function to compute the greeks ? 