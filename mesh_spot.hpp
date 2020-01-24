

#ifndef mesh_spot_HPP
#define mesh_spot_HPP

#include <vector>


namespace project
{
	class mesh
	{
	
	public: 
	
		mesh(const double& spot, const double& maturity, const double& volatility,const long& time_step,const long& steps, PayOff* _option);
		~mesh();
		PayOff* option;
		
		std::vector<double> Getvector_time() const;
		std::vector<double> Getvector_stock() const;
		double getdx() const;
		double getdt() const;
		double get_Spot() const;
		double init_cond(const double& x,const double& df = 0) const;
		std::vector<double> get_init_vector() const;
	
	
	private: 
	
		std::vector<double> m_vector_time;
		std::vector<double> m_vector_stock;
		double m_dx;
		double m_dt;
		double m_spot;
		long m_steps;
		long m_time_step;
		std::vector<double> m_init_vector;
		
	};
}

#endif