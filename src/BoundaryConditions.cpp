#include "BoundaryConditions.h"
using namespace std;



// Class constructor 
BoundaryConditions::BoundaryConditions(int id){ 
	this->id = id;

};

// Class destructor 
BoundaryConditions::~BoundaryConditions(void){
}; 

int BoundaryConditions::getId() const{
	return id;	
}

void BoundaryConditions::setId(int id){
	this->id = id;
}

string BoundaryConditions::getName(){
	return name;	
}

void BoundaryConditions::setName(string name){
	this->name = name;
}

string BoundaryConditions::getType(){
	return type;	
}

void BoundaryConditions::setType(string type){
	this->type = type;
}

int BoundaryConditions::getNoOfFaces() const{
	return noOfFaces;	
}

void BoundaryConditions::setNoOfFaces(int noOfFaces){
	this->noOfFaces = noOfFaces;
}

int BoundaryConditions::getStartFace() const{
	return startFace;	
}

void BoundaryConditions::setStartFace(int startFace){
	this->startFace = startFace;
}

void BoundaryConditions::setVariableType(int variableNo, string Type){
	this->variableType.setValue(variableNo, Type);
}

string BoundaryConditions::getVariableType(int variableNo){
	return this->variableType.getValue(variableNo);
}

TensorO1<string>* BoundaryConditions::getVariableTypeArray(){
	return &this->variableType;
};

void BoundaryConditions::setVariableValue(int variableNo, double Value){
	this->variableValue.setValue(variableNo, Value);
}

double BoundaryConditions::getVariableValue(int variableNo){
	return this->variableValue.getValue(variableNo);
}


TensorO1<double>* BoundaryConditions::getVariableValueArray(){
	return &this->variableValue;
};

TensorO1<string>* BoundaryConditions::getFunc2bcondTypeArray(){
	return &this->func2bcondType;
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Start of Boundary Conditions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


// A) Scalar Linear advection equation

// 1: fixedValue Boundary conditions
//____________________________________


void BoundaryConditions::fixedValueScalar(Face *face, int id, int rkStep){
	// fixedValue boundary conditions for scalarLinear solver
	TensorO1<double> fluxx(face->getNDOF());
	TensorO1<double> fluxy(face->getNDOF());
	TensorO1<double> fluxz(face->getNDOF());

	// Storing bcond[bcondi].variableValue[] in local variables

	double T = this->getVariableValue(0);
	TensorO1<double> vel(3);

	vel.setValue(0, this->getVariableValue(1));
	vel.setValue(1, this->getVariableValue(2));
	vel.setValue(2, this->getVariableValue(3));

	for (int DOF=0; DOF<face->getNDOF(); DOF++){
		fluxx.setValue(DOF, vel.getValue(0) * T);
		fluxy.setValue(DOF, vel.getValue(1) * T);
		fluxz.setValue(DOF, vel.getValue(2) * T);
	};

	// Normal flux at the face
	for (int DOF=0; DOF< face->getNDOF(); DOF++){
		double FLUX = face->getNormal()->getValue(0)*fluxx.getValue(DOF) + 
			      face->getNormal()->getValue(1)*fluxy.getValue(DOF) + 
			      face->getNormal()->getValue(2)*fluxz.getValue(DOF);
		face->setFlux(0,DOF, FLUX);
		// 0 in face.Flux[0][DOF]  stands for scalar variable
	};
};

// 2: Zero gradient Boundary conditions
//____________________________________

void BoundaryConditions::zeroGradientScalar(Face *face, int id, int rkStep){
	TensorO1<double> fluxx(face->getNDOF());
	TensorO1<double> fluxy(face->getNDOF());
	TensorO1<double> fluxz(face->getNDOF());
	TensorO1<double> vel(3);
	double T;


	for (int DOF=0; DOF<face->getNDOF(); DOF++){
		int DOFc = face->getOwnerDOFPoint(DOF);
		T = face->getOwnerCell()->getVariable(rkStep, 0, DOFc);
		vel.setValue(0, face->getOwnerCell()->getVariable(0, 1, DOFc));
		vel.setValue(1, face->getOwnerCell()->getVariable(0, 2, DOFc));
		vel.setValue(2, face->getOwnerCell()->getVariable(0, 3, DOFc));
		fluxx.setValue(DOF, vel.getValue(0) * T);
		fluxy.setValue(DOF, vel.getValue(1) * T);
		fluxz.setValue(DOF, vel.getValue(2) * T);
	};

		
	// Normal flux at the face
	for (int DOF=0; DOF< face->getNDOF(); DOF++){
		double FLUX = face->getNormal()->getValue(0)*fluxx.getValue(DOF) + 
			      face->getNormal()->getValue(1)*fluxy.getValue(DOF) + 
			      face->getNormal()->getValue(2)*fluxz.getValue(DOF);
		face->setFlux(0,DOF, FLUX);
		// 0 in face.Flux[0][DOF]  stands for scalar variable
	};
};








// B) Eulers equations of Gasdynamics

// 1: Fixed Value Boundary conditions
//____________________________________

void BoundaryConditions::fixedValueEuler(Face *face, int varId, int rkStep){
	// fixedValue boundary conditions for scalarLinear solver

	// Storing bcond[bcondi].variableValue[] in local variables

	for (int DOF=0; DOF<face->getNDOF(); DOF++){
		face->setPrimitiveVariable(varId,DOF, this->getVariableValue(varId));
	};
};


// 2: Zero gradient Boundary conditions
//____________________________________
void BoundaryConditions::zeroGradientEuler(Face *face, int varId, int rkStep){


	TensorO1<double> conservedVariableVector(5);
	TensorO1<double> primitiveVariableVector(5);
	for (int DOF=0; DOF<face->getNDOF(); DOF++){
		int DOFc = face->getOwnerDOFPoint(DOF); 
		for (int nVar=0; nVar<5; nVar++){
			conservedVariableVector.setValue(nVar, face->getOwnerCell()->getVariable(rkStep,nVar,DOFc));
		};
		getPrimitiveVariableVector(&conservedVariableVector, &primitiveVariableVector);

		face->setPrimitiveVariable(varId,DOF,primitiveVariableVector.getValue(varId));
	};
};


// 3: No Slip Boundary conditions
//____________________________________
void BoundaryConditions::noSlipEuler(Face *face, int varId, int rkStep){

	TensorO1<double> conservedVariableVector(5);
	TensorO1<double> primitiveVariableVector(5);
	TensorO1<double> fluxVectorX(5);
	TensorO1<double> fluxVectorY(5);
	TensorO1<double> fluxVectorZ(5);

	for (int DOF=0; DOF<face->getNDOF(); DOF++){
		
		//first fill the conserved variable vector
		int DOFc = face->getOwnerDOFPoint(DOF); 
		for (int nVar=0; nVar<5; nVar++){
			conservedVariableVector.setValue(nVar, face->getOwnerCell()->getVariable(rkStep,nVar,DOFc));
		};
		//get rho, u, v, w, p from the conserved variables
		getPrimitiveVariableVector(&conservedVariableVector, &primitiveVariableVector);
		// set velocity components to zero
		primitiveVariableVector.setValue(1, 0.0);
		primitiveVariableVector.setValue(2, 0.0);
		primitiveVariableVector.setValue(3, 0.0);

		face->setPrimitiveVariable(varId,DOF,primitiveVariableVector.getValue(varId));
	};
};

// 4: empty Boundary conditions
//____________________________________
void BoundaryConditions::emptyEuler(Face *face, int varId, int rkStep){
	// Normal flux at the face
	for (int DOF=0; DOF<face->getNDOF(); DOF++){
		face->setPrimitiveVariable(varId,DOF, 0.0);
	};
}; 










// C) Heat equation
void BoundaryConditions::zeroGradientHeat(Face *face, int id, int rkStep){
};
void BoundaryConditions::fixedValueHeat(Face *face, int id, int rkStep){
};


