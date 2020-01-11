#ifndef PAY_OFF_HPP
#define PAY_OFF_HPP

#include <algorithm> // Act on containers through iterators to apply modyfing/non_modifying operations

namespace projet
{
	//Our goal is to create an abstract class called PayOff, which will define an interface for all following payoff classes created by the user
	//Why doing this? This approach allows to encapsulate multiple different payoffs without modifying the existing ones

	//We use inheritance
	//Here is the base class
	class PayOff 
	{
		 public:
		  PayOff();
		  // Virtual destructor to avoid memory leaks when destroying the base and inherited classes
		  virtual ~PayOff() {}; 
		  // We turn the class into a functor (object we can call just like an object)
		  virtual double operator() (const double& S) const = 0;
	};
	//All the following class are bound to the Payoff class
	// first payoff class created for european call options
	class PayOffCall : public PayOff 
	{
		public:
		  PayOffCall(const double& K_);
		  virtual ~PayOffCall() {};

		  // Virtual function is now over-ridden (not pure-virtual anymore)
		  virtual double operator() (const double& S) const;
		  
		 private:
		  double K; // Variable for the Strike price

		 
	};

	class PayOffPut : public PayOff 
	{
		 private:
		  double K; // Strike

		 public:
		  PayOffPut(const double& K_);
		  virtual ~PayOffPut() {};
		  virtual double operator() (const double& S) const;
	};
}
#endif