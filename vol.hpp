//Partie hpp	
	//Classe de base, permet aussi le cas le plus simple ou la vol est constante
	
#ifndef VOL_HPP
#define VOL_HPP
#include <iostream>
#include <vector>
#include <algorithm>

#include "mesh_spot.hpp"

namespace project {
	class volatility 
	{
	
	public: 
		//Constructor
		volatility(const double& v,const mesh& grid); //Constructor, base class only takes initial value as input
		
		//Method to compute a vol (vol const donc les deux coeffs valent zéro)
		double compute_vol(const double& S = 0.,const double& t = 0.,const double& x = 0.,const double& y = 0.); //To compute a vol
		//Method to compute vector or vol (for one step of time)
		std::vector<double> vector_vol(const std::vector<double> v_spot,const double& t);
		
	private:
	
		std::vector<double> m_vol_const;

	
	protected:
		
		std::vector<double> m_vector_time;
		std::vector<double> m_vector_stock;
		std::size_t m_nb_time;
		std::size_t m_nb_spot;	
		double m_init_vol;
		mesh m_grid;
		
	};
	 //Classe ou on va dire que la volatility est fonction du temps et du spot
	//On pourra créer ainsi d'autres classes avec différents comportement pour la vol (e.g. seulement une fonction du temps)
	class vol_surface : public volatility
	{
	public:
	
		//Niveau paramètre, la classe prends en plus les deux coeff pour définit la fonction
		vol_surface(const double& v, mesh grid, const double& coeff_tps, const double& coeff_spot);
		
		//On reprends les mêmes méthodes (du coup la les coeffs ne sont plus constant zéro comme au dessus)
		double compute_vol(const double& S,const double& t,const double& x,const double& y); 
		std::vector<double> vector_vol(const std::vector<double> v_spot,const double& t);

		
	private:
	
		double m_coeff_time;
		double m_coeff_spot;
		std::vector<double> m_vol_matrix;
		
	
	};
}
#endif