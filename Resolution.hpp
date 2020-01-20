#ifndef Resolution_hpp
#define Resolution_hpp

#include <vector>
#include <cmath>
#include <iostream>
//Allows to create and manipulate xarrays
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xeval.hpp"
//Allows to perform linear algebra operations on xarrays
#include "xtensor-blas/xlinalg.hpp"

class Solve
    {
        
    public:

        Solve(mesh _grid, Parameters _param, bound_conditions _conditions);
		
		~Solve();
		
		
        
        //double dt; //to get the time step
        //double dx; //to get the stock value step
        //int nb_step_time; //to get the number of time steps on which we iterate
        //double v; //to get the volatility
        //double r; //to get the rate
        //double theta; //to get the theta
        //int nb_step_spot; //to get the size of the matrix A and B (N-1)


        //Iterates on each time step of the grid to solve the equation. Beginning from time T-1 and solving for time 1

    };
#endif /* Resolution_hpp */
