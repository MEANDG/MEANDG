#include<iostream>
#include"InternalFluxSolver.h"

using namespace std;

void scalarTransportFlux(Cell *cell, int rkstep){

	// analytical flux at DOF locations
	double vel[3];
	double T;

	int nDOF = cell->getNDOF();
	// velocity vector is not getting updated. So 0th instant value is copied
	for (int DOF=0; DOF< nDOF; DOF++){
		T = cell->getVariable(rkstep,0,DOF);
		vel[0] = cell->getVariable(0,1,DOF);
		vel[1] = cell->getVariable(0,2,DOF);
		vel[2] = cell->getVariable(0,3,DOF);

		cell->getFluxVectorx()->setValue(0, DOF, vel[0] * T); // 0 indicates scalar equation
		cell->getFluxVectory()->setValue(0, DOF, vel[1] * T); // 0 indicates scalar equation
		cell->getFluxVectorz()->setValue(0, DOF, vel[2] * T); // 0 indicates scalar equation
	};
};


void BurgerFlux(Cell *cell, int rkstep){
	// analytical flux at DOF locations
	double vel[3];
	double T;

	int nDOF = cell->getNDOF();
	for (int DOF=0; DOF< nDOF; DOF++){
			   // Burger equation //
		T = cell->getVariable(rkstep,0,DOF);
		cell->getFluxVectorx()->setValue(0,DOF, sqr(T)/2.0);
		cell->getFluxVectory()->setValue(0,DOF, 0.0);
		cell->getFluxVectorz()->setValue(0,DOF, 0.0);
	};
};


void HeatEquationFlux(Cell *cell, int rkstep){

};


void EulerFlux(Cell *cell, int rkstep){
	// analytical flux at DOF locations
	// note that although the initial conditions are given as a primitive variable vector, when we distribute IC on DOF locations,
	// we first create a conserved variable vector and then copy this vector for cell DOF locations. 
	// thus DG::variable has primitive variables but Cell::variable has conserved variables. This has to be kept in mind while 
	// writing the output

	TensorO1<double> conservedVariableVector(5);
	TensorO1<double> FluxVectorx(5);
	TensorO1<double> FluxVectory(5);
	TensorO1<double> FluxVectorz(5);

	int nDOF = cell->getNDOF();
	// velocity vector is not getting updated. So 0th instant value is copied
	for (int DOF=0; DOF< nDOF; DOF++){

		// copy values of the conserved variable Vector
		for (int varNo=0; varNo<5; varNo++){
			conservedVariableVector.setValue(varNo,  cell->getVariable(rkstep,varNo,DOF));
			// reset Flux to zero
			FluxVectorx.setValue(varNo, 0.0);
			FluxVectory.setValue(varNo, 0.0);
			FluxVectorz.setValue(varNo, 0.0);
		};
		

		// Get the actual Flux at DOF locations 
		getEulerFluxFromConservedVariables(&conservedVariableVector, &FluxVectorx, &FluxVectory, &FluxVectorz ); // refer src/gasdynamics.cpp

		// copy  this for cell->fluxVectors 
		for (int varNo=0; varNo<5; varNo++){
			cell->getFluxVectorx()->setValue(varNo,DOF, FluxVectorx.getValue(varNo));
			cell->getFluxVectory()->setValue(varNo,DOF, FluxVectory.getValue(varNo));
			cell->getFluxVectorz()->setValue(varNo,DOF, FluxVectorz.getValue(varNo));
		};
	};
	//done
};

void NavierStokesFlux(Cell *cell, int rkstep){
};


