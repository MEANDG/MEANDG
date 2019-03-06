#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H
#include "sysInclude.h"
#include "mathFunctions.h"
#include "Point.h"
#include "Face.h"
#include "Cell.h"
#include "Gasdynamics.h"

class BoundaryConditions{
	
public:
	/// Constructor function
	BoundaryConditions(int id = 0);				
	
	/// Destructor function
	~BoundaryConditions(); 

	/// Get method for boundary id
	int getId() const;

	/// Set method for boundary id
	void setId(int id);

	/// Get the boundary name
	string getName();

	/// Set the boundary name
	void setName(string name);
	
	/// Get the boundary type
	string getType();

	/// Set the boundary type
	void setType(string type);

	/// Get the number of faces
	int getNoOfFaces() const;

	/// Set the number of faces
	void setNoOfFaces(int noOfFaces);
	
	/// Get the starting face
	int getStartFace() const;

	/// Set the starting face
	void setStartFace(int startFace);

	/// Set variable type for all the primitive variables
	void setVariableType(int variableNo, string Type);

	/// Set variable type for all the primitive variables
	string getVariableType(int variableNo);
	
	/// Get a pointer to the variableType array
	TensorO1<string>* getVariableTypeArray();

	/// Set variable value for all the primitive variables
	void setVariableValue(int variableNo, double Value);

	/// Set variable type for all the primitive variables
	double getVariableValue(int variableNo);

	/// Get a pointer to the variableValue array
	TensorO1<double>* getVariableValueArray();

	/// Get a pointer to the func2bcondType array
	TensorO1<string>* getFunc2bcondTypeArray();

	//void (BoundaryConditions::*func2bcond)(Face * face, int rkstep);
	/// Function pointer pointing at the appropriate boundary condition
	vector <void (BoundaryConditions::*) (Face *, int, int)> func2bcond;				

	/*__________________________________________________________________________________________*/
	//  1. Boundary conditions for scalar problem

	/// fixedValue boundary condition for the scalar problem
	void fixedValueScalar(Face *face, int id, int rkStep);

	/// zeroGradient boundary conditions for the scalar problem 
	void zeroGradientScalar(Face *face, int id, int rkStep);

	/*__________________________________________________________________________________________*/
	//  2. Boundary conditions for Euler's equations
	
	/// zeroGradient boundary condition for Euler's equation
	void zeroGradientEuler(Face *face, int id, int rkStep);

	/// fixedValue boundary condition for Euler's equations
	void fixedValueEuler(Face *face, int id, int rkStep);

	/// noSlip (solid wall) boundary condition for Euler's equations
	void noSlipEuler(Face *face, int id, int rkStep);

	/// empty (zeroGrad) boundary condtions for Euler's equations
	void emptyEuler(Face *face, int id, int rkStep);

	/*__________________________________________________________________________________________*/
	//  3. Boundary conditions for Heat equations
	
	// euler equations boundary conditions
	void zeroGradientHeat(Face *face, int id, int rkStep);
	void fixedValueHeat(Face *face, int id, int rkStep);

	/*__________________________________________________________________________________________*/
private:
	/// OpenFOAM identifier of the boundary
	int id;

	/// Boundary name
	string name;						

	/// Boundary type 
	string type;						

	/// Number of faces for the boundary
	int noOfFaces;						

	/// Starting face for the boundary 
	int startFace;						

	/// Variable type for each of the primitive variables for this boundary condition (1 type for scalar)
	TensorO1<string> variableType;	

	/// Value of the primitive variable for this boundary condition
	TensorO1<double> variableValue;

	/// Function type assigned for each entry in func2bcond (i.e. fixedValueEuler if system type is Euler and boundary condition type is fixedValue. Assigned for all the variables for all the Boundary conditions)
	TensorO1<string> func2bcondType;

};

#endif
