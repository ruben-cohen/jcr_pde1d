


#ifndef BOUND_CONDITIONS_HPP
#define 
#include <vector>

// we create a boundaries conditions class that is purely virtual and 2 other classes neumann and derichtlet


namespace project{
	
// boundaries class
	
	class bound_conditions 
	{
		public:
		
			bound_conditions(mesh _grid);
			
			//std::vector<double> cond() const;
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 //std::vector<double>  cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
		
		protected:
		
			mesh m_grille;
			//Parameters m_param;
	
		//private:
		
		//std::vector<std::vector<double>> Matrix_conditions;
			
	};
	class Derichtlet: public bound_conditions 
	{
		
		public:
		
			Derichtlet(const mesh& m_grid, const std::vector<double>& rate);
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);
			 //std::vector<double> cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 std::vector<double> get_cond() const;

			std::vector<double> matrix_derichtlet;
			
		//std::vector<double> K_neuman; 
		
	};
	
	
	class Neumann: public bound_conditions 
	{
		
		public:
		
			Neumann(mesh m_grid, double theta, std::vector<double> sigma, std::vector<double> rate,std::vector<double>& const_vector);
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);
			 //std::vector<double> cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 std::vector<double> get_cond() const;

			std::vector<double> matrix_neumann;
			
		//std::vector<double> K_neuman; 
		
	};
		
}


#endif 