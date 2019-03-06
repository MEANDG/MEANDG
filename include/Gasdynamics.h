#ifndef GASDYNAMICS_H
#define GASDYNAMICS_H


#include "sysInclude.h"
#include "mathFunctions.h"

void getConservedVariableVector(TensorO1<double> *primitiveVariableVector, TensorO1<double> *conservedVariableVector);
void getPrimitiveVariableVector(TensorO1<double> *conservedVariableVector, TensorO1<double> *primitiveVariableVector);
void getEulerFluxFromConservedVariables(TensorO1<double> *conservedVariableVector, TensorO1<double> *fluxVectorX, TensorO1<double> *fluxVectorY, TensorO1<double> *fluxVectorZ);
void getEulerFluxFromPrimitiveVariables(TensorO1<double> *primitiveVariableVector, TensorO1<double> *fluxVectorX, TensorO1<double> *fluxVectorY, TensorO1<double> *fluxVectorZ);

void findRoeAverages(TensorO1<double> *leftVector, TensorO1<double> *rightVector, TensorO1<double> *avgVector);

void findEigenValuesForEuler(TensorO1<double> *conservedVariableVector, TensorO1<double> *normal, TensorO1<double> *Lambda);
void findEigenVectorsForEuler(TensorO1<double> *conservedVariableVector, TensorO1<double> *normal, TensorO1<double> *tangent1, TensorO1<double> *tangent2, Matrix<double> *Rn, Matrix<double> *Ln);

double expp(double x, double y); 
void getAnalyticalShockTube(TensorO2<double> *SodAnalytical);

void getNSFluxFromConservedVariables(TensorO1<double> *conservedVariableVector, TensorO1<double> *fluxVectorX, TensorO1<double> *fluxVectorY, TensorO1<double> *fluxVectorZ, double nu);


//**********************************************************************************************************************//
// TESTS, Tests
//**********************************************************************************************************************//

namespace Test{
	void TestGasdynamicsFunctions();
};

#endif
