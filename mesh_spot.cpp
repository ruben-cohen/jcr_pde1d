
#include 'mesh_spot.hpp'
#include <vector>
#include <cmath>



namespace project{
	
	
	mesh_spot::grid_mesh(double spot, double maturity, double volatility, double time_step, size_t steps){
		: dx(dx_spot), dt(time_steps)
		
		
		double stock_init = log(spot);
		double high_bound = stock_init + 5*volatility*sqrt(maturity);
		double low_bound = stock_init - 5*volatility*sqrt(maturity);
		
		double dx_spot = (high_bound - low_bound)/steps;
		
		std::vector<double> vector_stock(steps);
		
		for (std::size_t i = 0; i < steps ; ++i){
			
			vector_stock[i] = stock_init + (i - steps )* dx_spot;

		}
		
		std::vector<double> vector_time(time_step);
		
		for (std::size_t j = 0; j < maturity ; ++j){
			
			vector_time[i] = j*time_steps;

		}
	};
	mesh_spot::~grid_mesh() {

	} //destructor of the grid 
	
	std::vector<double> mesh_spot::Getvector_time()const{
		
			return vector_time;
	} //useful to get the vector of time from the mesh 
			
	std::vector<double> mesh_spot::Getvector_stock() const{
		
		return vector_stock;
	
	} //useful to get the vector of stock path from the mesh 
			
	double mesh_spot::getdx() const{
		
		return dx_spot;
	
	} //we will need to get the dx for CN algo 
	
	double mesh_spot::getdt() const{
		
		return time_steps;
	
	} // we will need the dt for CN algo 
	
		
		
		
		
		
		
		
		
	}
	
	
	
	
	
}