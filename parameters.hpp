


#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>

// this class is only the initialization of the three parameters needed to create the matrix 

namespace project{
	
		class Parameters {
	public:
		Parameters(double vol, double rate, double theta);
		double Get_Vol() const;
		double Get_Rate() const;
		double Get_Theta() const;
		~Parameters();

	private:
		
		double pa_vol;
		double pa_Rate;
		double pa_Theta;
	};
		
}


#endif 