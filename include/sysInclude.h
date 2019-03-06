/** @file */
/// Includes all the systemwide header files 
#include<iostream>
#include<cstdlib>
#include<string>
#include<fstream>
#include<sstream>
#include<assert.h>
#include<iomanip>
#include<sys/stat.h>
#include<cmath>
#include<math.h>
#include<time.h>
#include<vector>
#include <omp.h>

//Default namespace is std for this program
using namespace std;

#ifndef CONSTANTS
#define CONSTANTS
const double PI = 3.14159265;
const double GAMMA = 1.4; 	// specific heat ratio

/// Tolerances
const double EPS = 2.220446049250313e-16;		//machine accuracy

const double SMALL = pow(10,-3);
const double MICRO = pow(10,-6);	
const double NANO = pow(10,-9);
const double PICO = pow(10,-12);
const double FEMTO = pow(10,-15);
const double SMALL10 = pow(10,-10);

const double LARGE= pow(10,3);
const double MEGA = pow(10,6);
const double GIGA = pow(10,9);
const double TERA = pow(10,12);
const double PETA = pow(10,15);
const double LARGE10 = pow(10,10);
#endif

/// Define dTypes
typedef unsigned int Index;
typedef unsigned int unt;

// Include all the other supporting files

// 0. Functions which do not depend on anything 
#include "cellTypes.h"
#include "faceTypes.h"
#include "Eigen/Dense" 	//include Eigen library for matrix operations
#include "intFlag.h" 
#include "Solvers.h"
#include "Systems.h"

// 1. Functions which do not depend on other MEANDG functions but may depend on sysInclude.h functions
#include "displayFunctions.h"
#include "otherFunctions.h"
#include "randomGenerator.h"
#include "Tensor.h"   



