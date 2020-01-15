

#include 'bound_conditions.hpp'
#include<math>

namespace project {
	
	
	std::vector<double> operator() Neumann::neumann(mesh grid, parameters param, Payoff* option,std::vector<double> K_neuman){
	
	double dx = grid.getdx(); //need the stock step 
	double sigma = param.get_Vol(); //to get the volatility 
	double rate = param.get_Rate(); //the get the rate 
	double theta = param.get_Theta(); //to get the theta 
	double size_vec = grid.Getvector_time().size();
	double maturity = grid.Getvector_time().back();
	double S0 = grid.get_Spot();
	
	double K1 = K_neuman[0];
	double K2 = K_neumann[1];
	double K3 = K_neumann[2];
	double K4 = K_neumann[3];
	
	PDE _payoff(option,dx,dt,grid.Getvector_time(),grid.Getvector_stock()); //create the PDE object from the option 
	std::vector<double>  init_cond = _Payoff.init_cond(S0); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
	
	double f_0_T = init_cond[0]; //first element of the vector is the payoff at min S and maturity 
	double f_N_T = init_cond.back(); //last element is the payoff at max S and maturity 
	
	std::vector<double> upper_conditions(size_vec); 
	std::vector<double> lower_conditions(size_vec); 
	
	std::fill (upper_conditions.begin(),upper_conditions.end()-1,0);   // we fill the vector with 0 at time 1 to T-1 
	upper_conditions.back() = f_0_T*exp(-maturity*rate);
	std::fill (lower_conditions.begin(),lower_conditions.end()-1,0);
	lower_conditions.back() = f_N_T*exp(-maturity*rate);
	
	//double beta_left =  dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
	//double beta_right = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
	
	//double alpha_left = dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
	//double alpha = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
	
	double coef_ =(1 - dt*(1-theta)*rate)/(1+ dt*theta*rate);
	double coef_K1_K2 = -dt*((-sigma**2)*K1)/2 + ((sigma**2)/2 - rate)*K2))/(1+ dt*theta*rate);
	double coef_K3_K4 = -dt*((-sigma**2)*K3)/2 + ((sigma**2)/2 - rate)*K4))/(1+ dt*theta*rate);
	
	for (auto it = upper_conditions.rbegin(); it =! upper_conditions.rend(); it++){
	//rbegin() reverse iterator to fill the vector from the end to the beginning
		
		upper_conditions[it] = coef_*upper_conditions[it-1] + coef_K1__K2;
		
		lower_conditions[it] = coef_*lower_conditions[it-1] + coef_K3__K4;
	}

		
	};
	
	 
	std::vector<double> operator() Derichtlet::derichtlet(mesh grid, parameters param, Payoff* option,std::vector<double> K_neuman ={0,0,0,0}){
		
	 
	double dx = grid.getdx(); //need the stock step 
	double sigma = param.get_Vol(); //to get the volatility 
	double rate = param.get_Rate(); //the get the rate 
	double theta = param.get_Theta(); //to get the theta 
	double size_vec = grid.Getvector_time().size();
	double maturity = grid.Getvector_time().back();
	double S0 = grid.get_Spot();
	
	
	PDE _payoff(option,dx,dt,grid.Getvector_time(),grid.Getvector_stock()); //create the PDE object from the option 
	std::vector<double>  init_cond = _Payoff.init_cond(S0); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
	
	double f_0_T = init_cond[0]*exp(-maturity*rate); //first element of the vector is the payoff at min S and maturity 
	double f_N_T = init_cond.back()*exp(-maturity*rate); //last element is the payoff at max S and maturity 
	
	std::vector<double> upper_conditions(size_vec); 
	std::vector<double> lower_conditions(size_vec); 
	
	std::fill (upper_conditions.begin(),upper_conditions.end(),f_0_T);   // we fill the vector with the terminal condition as Vn f_n is a constant (l or h)  
	std::fill (lower_conditions.begin(),lower_conditions.end(),f_N_T);
	
	//double beta_left =  dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
	//double beta_right = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (rate)/(2*dx);
	
	//double alpha_left = dt*theta*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
	//double alpha = -dt*(1-theta)*(-sigma**2)/(2*dx**2) + (sigma**2)/(4*dx**2) + (-rate)/(2*dx);
	
	double coef_ =1;
	double coef_K1_K2 = 0;
	double coef_K3_K4 = 0;
	
	for (auto it = upper_conditions.rbegin(); it =! upper_conditions.rend(); it++){
	//rbegin() reverse iterator to fill the vector from the end to the beginning
		
		upper_conditions[it] = coef_*upper_conditions[it-1] + coef_K1__K2;
		
		lower_conditions[it] = coef_*lower_conditions[it-1] + coef_K3__K4;
	}

		
	};
	
	std::vector<double>  bound_conditions::boundaries_compute(mesh grid, parameters param, Payoff* option, bound_conditions* bound_func, std::vector<double> K_neuman ={0,0,0,0}){
		
	return (*bound_func)(mesh_spot grid, parameters param, Payoff* option, std::vector<double> K_neuman ={0,0,0,0});
		
	}
	
	std::vector<double> bound_conditions::get_upper_border(){
		
		return upper_conditions;
	};
	
	std::vector<double> bound_conditions::get_lower_border(){
		
		return lower_conditions;
	};
	
	
	bound_conditions::~bound_conditions(){}; //destructor of the boundaries_conditions 
	
	
	
	
	
	
	
	
	
	
	
	
	
}