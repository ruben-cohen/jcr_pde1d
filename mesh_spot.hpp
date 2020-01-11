

#ifndef mesh_spot_HPP
#define mesh_spot_HPP

#include <vector>


namespace project
{
	class mesh{
	
	public: 
	grid_mesh(double spot, double maturity, double volatility, double time_step, size_t steps);
	~grid_mesh();
	std::vector<double> Getvector_time() const;
	std::vector<double> Getvector_stock() const;
	double getdx() const;
	double getdt() const;
	double get_Spot const;
	
	
	private: 
	
	std::vector<double> vector_time;
	std::vector<double> vector_stock;
	double dx;
	double dt;
	double spot;
	
	
	}  
}

#endif