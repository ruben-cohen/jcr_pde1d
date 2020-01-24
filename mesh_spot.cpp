
#include 'mesh_spot.hpp'
#include <vector>
#include <cmath>



namespace project{
	
	
	//This class builds the mesh, giving dt, dS, and the axes of the mesh as output
	mesh::mesh(const double& spot, const double& maturity,const double& volatility, const long& time_step, const long& steps,PayOff* _option)
	:
	 option(_option),
	 m_dt(maturity/time_step),
	 m_spot(spot),
	 m_steps(steps),
	 m_time_step(time_step)
	{
		//The mesh is centered around log(S0)
		double S0 = log(m_spot);
		double high_bound = S0 + 5*volatility*sqrt(maturity);
		double low_bound = S0 - 5*volatility*sqrt(maturity);
		
		//dS
		m_dx = (high_bound - low_bound)/steps;
		
		//Axe of spot price
		for (long i = 0; i < m_steps+1 ; ++i)
		{
			
			m_vector_stock.push_back(high_bound + (i - m_steps)*m_dx);
		}
		//Axe of time
		for (std::size_t j = 0; j < m_time_step+1 ; ++j)
		{
			m_vector_time.push_back(j*m_dt);	
		}
		
		size_t nb_step = m_vector_stock.size();
		//std::vector<double> m_init_vector(nb_step);
		
		for (std::size_t i = 0; i < nb_step ; ++i)
		{
			m_init_vector.push_back(init_cond(exp(m_vector_stock[i])));
		}
	};

		
	//Multiple Get methods to return private variables of the class


	//useful to get the vector of time from the mesh 
	std::vector<double> mesh::Getvector_time()const
	{
		return m_vector_time;
	}; 
	//useful to get the vector of stock path from the mesh 
	std::vector<double> mesh::Getvector_stock() const
	{
		return m_vector_stock;
	}; 
	double mesh::init_cond(const double& x,const double& df) const 
	{
	  return option->operator()(x,df);
	};
	
	//we will need to get the dx for CN algo	
	double mesh::getdx() const
	{
		return m_dx;
	}; 
	std::vector<double> mesh::get_init_vector() const //const forbid to modify the state of my object
    {
        return m_init_vector;
    };
	
	// we will need the dt for CN algo 
	double mesh::getdt() const
	{
		return m_dt;
	}; 
	
	double mesh::get_Spot() const
	{
		return m_spot;
	};
	
	//Destructor
	mesh::~mesh() {};
}
	
