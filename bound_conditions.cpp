

#include 'bound_conditions.hpp'
#include<math>

namespace project {
	
	
	std::vector<double> operator() Neumann::neumann(mesh_spot grid, parameters param, Payoff *function_payoff,std::vector<double> K_neuman){
	
	double dx = grid.getdx(); //need the stock step 
	double sigma = param.get_Vol(); //to get the volatility 
	double rate = param.get_Rate(); //the get the rate 
	double theta = param.get_Theta(); //to get the theta 
	double size_vec = grid.Getvector_time().size();
	double maturity = grid.Getvector_time.back();
	double  spot = grid.get_Spot();
	
	double K1 = K_neuman[0];
	double K2 = K_neumann[1];
	double K3 = K_neumann[2];
	double K4 = K_neumann[3];
	
	
	std::vector<double> right_conditions(size_vec); 
	std::vector<double> left_conditions(size_vec); 
	
	std::fill (left_conditions.begin()+1,left_conditions.end()-1,0);   // we fill the vector with 0 at time 1 to T-1 
	std::fill (right_conditions.begin()+1,right_conditions.end()-1,0);
	
	double beta_left =  dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
	double beta_right = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
	
	double alpha_left = dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
	double alpha_right = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
	
	double coef_left_0 =   (1 - dt*(1-theta)*rate)/(1+ dt*theta*rate);
	double coef_K1_K2 = -dt*((-sigma**2)*K1)/2 + ((sigma**2)/2 - rate)*K2))/(1+ dt*theta*rate);
	double coef_K3_K4 = -dt*((-sigma**2)*K3)/2 + ((sigma**2)/2 - rate)*K4))/(1+ dt*theta*rate);
	
	// sum cumulative mais dépend du time step ? :/  
		
		
	double right_begin = beta_right*0; 
	double right_end = alpha_right*function_payoff(spot)* exp(-maturity*rate); 
	
	left_conditions.front() = left_begin;
	left_conditions.back() = left_end;
	
	right_conditions.front() = right_begin;
	right_conditions().back() = right_end;

	
		
	};
	
	 
	std::vector<double> operator() Derichtlet::derichtlet(mesh_spot grid, parameters param, Payoff *function_payoff,std::vector<double> K_neuman ={0,0,0,0}){
		
	 
	double dx = grid.getdx(); //need the stock step 
	double sigma = param.get_Vol(); //to get the volatility 
	double rate = param.get_Rate(); //the get the rate 
	double theta = param.get_Theta(); //to get the theta 
	double size_vec = grid.Getvector_time().size();
	double maturity = grid.Getvector_time.back();
	double  spot = grid.get_Spot();
	
	
	std::vector<double> right_conditions(size_vec); 
	std::vector<double> left_conditions(size_vec); 
	
	std::fill (left_conditions.begin()+1,left_conditions.end()-1,0);   // we fill the vector with 0 at time 1 to T-1 
	std::fill (right_conditions.begin()+1,right_conditions.end()-1,0);
	
	double beta_left =  dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
	double beta_right = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
	
	double alpha_left = dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
	double alpha_right = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
	
	double left_begin = beta_left*0; 
	double left_end = alpha_left*function_payoff(spot)* exp(-maturity*rate); 
		
		
	double right_begin = beta_right*0; 
	double right_end = alpha_right*function_payoff(spot)* exp(-maturity*rate); 
	
	left_conditions.front() = left_begin;
	left_conditions.back() = left_end;
	 
	right_conditions.front() = right_begin;
	right_conditions().back() = right_end;
	// sinon il faut le faire en récursif, mais je sais pas trop comment on récupère le bon time step ? 	
	};
	
	std::vector<double>  bound_conditions::boundaries_compute(mesh_spot grid, parameters param, Payoff* function_payoff, bound_conditions* bound_func, std::vector<double> K_neuman ={0,0,0,0}){
		
	return (*bound_func)(mesh_spot grid, parameters param, Payoff *function_payoff, std::vector<double> K_neuman ={0,0,0,0});
		
	}
	
	std::vector<double> bound_conditions::get_right_border(){
		
		return right_conditions;
	};
	
	std::vector<double> bound_conditions::get_left_border(){
		
		return left_conditions;
	};
	
	
	bound_conditions::~bound_conditions(){}; //destructor of the boundaries_conditions 
	
	
	
	
	
	
	
	
	
	
	
	
	
}