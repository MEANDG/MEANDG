#ifndef OTHERFUNCTIONS_H
#define OTHERFUNCTIONS_H

#include "sysInclude.h"

using namespace std;

	/// Converts a given datatype to a string;
	template<typename var>
	string varToString (var input){
		stringstream ss;
		ss << input;
		string output;
		output= ss.str();
		return output;
	};

	int stringToInt(string input);
	float stringToFloat(string input);
	double stringToDouble(string input);
	// Factorial function
	int fact(int k);


	/// Returns absolute value 
	template<typename T>
	T abs(T a)
	{
		if (a > 0)
		{
			return a;
		}
		else
		{
			return -a;
		}
	};


	/// Returns square of a number 
	template<typename T>
	T sqr(T a)
	{
		T b;
		b = a*a;
		return b;
	};

	
	bool fileExists(const std::string& filename);


	
	void generatePlottingScript(string name, int order);
	void generatePlottingScriptAssignment1(string name);
	void generatePlottingScriptAssignment2(string name);


	void generatePlottingScript2D(string name, int order);
#endif


