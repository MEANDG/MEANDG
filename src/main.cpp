/** @file */
#include "sysInclude.h"
#include "Tests.h"
using namespace std;


int main(){

	srand(int(time(0)));

	//runAllTests();
	Test::TestEulerEquations1();

	return 0;
}
