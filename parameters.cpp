
#include 'parameters.hpp'

// this class is only the initialization of the three parameters needed to create the matrix 
namespace project
{
	Parameters::Parameters(double vol, double rate, double theta)
		: pa_vol(vol), pa_Rate(rate), pa_Theta(theta)
	{

	}
	double Parameters::Get_Vol() const{

		return pa_vol;
	}
	double Parameters::GetRate() const{

		return pa_rate;
	}
	double Parameters::Get_Theta() const{

		return pa_theta;
	}
	Parameters::~Parameters() {

	}
}