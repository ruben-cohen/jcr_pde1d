

#include 'bound_conditions.hpp'
#include<math>
#include 'payoff.hpp'
#include 'mesh_spot.hpp"


//Boundaries	
	
//Boundaries	
	
	bound_conditions::bound_conditions(mesh _grid)
	:
	 m_grille(_grid)
	 //,m_param(_param)
	 {};
	 
	
	Derichtlet::Derichtlet(const mesh& m_grid,const std::vector<double>& rate)
	:
	 bound_conditions(m_grid)
	{
		//From mesh object
		double dt = m_grid.getdt();
		double dx = m_grid.getdx(); //need the stock step 
		size_t size_vec = m_grid.Getvector_time().size();
		double T = m_grid.Getvector_time().back();
		double maturity = m_grid.Getvector_time().back();
		double S0 = m_grid.get_Spot();
		size_t size_spot = m_grid.Getvector_stock().size();
		
/* 		if(size_spot != rate.size()){
			
			std::cout<< "size of sigma is" << rate.size() << "size of spot is" << size_spot << std::endl;
		}
		else{ std::cout<< "okay for Derichtlet" << std::endl;}; */
		
		//From parametre object
		//double sigma = m_param.Get_Vol(); //to get the volatility 
		//double rate = m_param.Get_Rate(); //the get the rate 
		//double theta = m_param.Get_Theta(); //to get the theta 
		//double r = m_param.Get_Rate();
		
		//From PDE object
		double spot_max = m_grid.Getvector_stock().back();
		double spot_min = m_grid.Getvector_stock()[0];
		//double f_0_T = _init_cond[0];
		//double f_N_T = _init_cond.back();
		
		matrix_derichtlet.resize(size_vec*2); 
		//std::vector<double> lower_conditions(size_vec);

		double df = 0.;
		
		for(size_t i = 0; i < size_vec; ++i)
        {
			df = exp(-rate[i]*(T-dt*i));
			matrix_derichtlet[i] = m_grid.init_cond(exp(spot_min),df); // for here
			matrix_derichtlet[size_vec+i] = m_grid.init_cond(exp(exp(spot_max),df);;
        }
	};
	
	std::vector<double> Derichtlet::get_cond() const
	{
		return matrix_derichtlet;
	};
	
	Neumann::Neumann(mesh m_grid, double theta, std::vector<double> sigma, std::vector<double> rate,std::vector<double>& const_vector)
	:bound_conditions(m_grille)
	{
		//From mesh object
		double dt = m_grid.getdt();
		double dx = m_grid.getdx(); //need the stock step 
		size_t size_vec = m_grid.Getvector_time().size();
		double maturity = m_grid.Getvector_time().back();
		double S0 = m_grid.get_Spot();
		size_t size_spot = m_grid.Getvector_stock().size();		
		
		//From parametre object
		//double sigma = m_param.Get_Vol(); //to get the volatility 
		//double rate = m_param.Get_Rate(); //the get the rate 
		//double theta = m_param.Get_Theta(); //to get the theta 
		//double r = m_param.Get_Rate();
		
/* 		if(size_spot != sigma.size()){
			
			std::cout<< "size of sigma is" << sigma.size() << "size of spot is" << size_spot << std::endl;
		}
		else{ std::cout<< "okay for neumann" << std::endl;};
		
		if(size_spot != rate.size()){
			
			std::cout<< "size of rate is" << rate.size() << "size of spot is" << size_spot << std::endl;
		}
		else{ std::cout<< "okay for neumann" << std::endl;}; */
		
		//From PDE object
		std::vector<double>  _init_cond = m_grid.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		double f_0_T = _init_cond[0];
		double f_N_T = _init_cond.back();
		
	    double K1 = const_vector[0];
	    double K2 = const_vector[1];
	    double K3 = const_vector[2];
	    double K4 = const_vector[3];
		   
		double coef_; 
		double coef_K1_K2; 
		double coef_K3_K4; 
		
		matrix_neumann.resize(size_vec);
		matrix_neumann.back() =  f_0_T*exp(-maturity*rate[0]);
		matrix_neumann.resize(size_vec*2);
		matrix_neumann.back() = f_N_T*exp(-maturity*rate.back());
		
		for (size_t it = size_vec-1; it > 0; it--)
		{
		      //reverse iterator to fill the vector from the end to the beginning
			  
			  //std::cout << it << std::endl;
			  
			  size_t tt = size_vec + it -1;
			  
			  //std::cout << tt << std::endl;
			  
			  		coef_ =(1 - dt*(1-theta)*rate[it])/(1+ dt*theta*rate[it]);
					coef_K1_K2 = -dt*((-std::pow(sigma[it],2)*K1)/2 + (std::pow(sigma[it],2)/2 - rate[it])*K2)/(1+ dt*theta*rate[it]);
					coef_K3_K4 = -dt*((-std::pow(sigma[it],2)*K3)/2 + (std::pow(sigma[it],2)/2 - rate[it])*K4)/(1+ dt*theta*rate[it]);
					
					
			  matrix_neumann[it]  = (coef_*matrix_neumann[it-1] + coef_K1_K2)*exp(-rate[it]*dt*it);
			
			  matrix_neumann[tt] = (coef_*matrix_neumann[tt-1] + coef_K3_K4)*exp(-rate[it]*dt*it);
		 }
		 
	};
	
    std::vector<double> Neumann::get_cond() const
	{
		return matrix_neumann;
	};
/////////////////////////////////////////////////
	
	
	