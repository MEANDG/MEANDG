#ifndef GETRANDOM
#define GETRANDOM

#include "sysInclude.h"

template<typename Type>
Type getRandom(Type Min, Type Max)
{

	double dMin = Min * 1.0;
	double dMax = Max * 1.0;
	double randomValue = (double)rand() / RAND_MAX;
	double value = Min + randomValue * (Max - Min);
	Type output = value;
	return output;
}
#endif

