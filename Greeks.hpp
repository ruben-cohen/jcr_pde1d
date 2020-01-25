#ifndef Greeks_hpp
#define Greeks_hpp
#include <iostream>
#include <vector>
#include <algorithm>
#include "mesh_spot.hpp"
#include "bound_conditions.hpp"

class Greeks
{
    
public:

    Greeks(mesh grid, solver, Neumann, Derichtlet);
    
    ~Greeks();
    
    
};




#endif /* Greeks_hpp */
