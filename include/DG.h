#ifndef DG_H
#define DG_H
#include "sysInclude.h"
#include "Point.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryConditions.h"
#include "GeometryIO.h"
#include "RefCell.h"
#include "RefFace.h"
#include "InternalFluxSolver.h"
#include "RiemannSolver.h"

class DG{
	
public:
// Global variables

// 0. Flags
	// Flags for memory management
	// Refer the destructor function of this class
	// Flag returns 'true' if memory is allocated using 'new' method. Then delete is required to free the memory.
	bool VariableNameFlag;
	bool pointsFlag;
	bool facesFlag;
	bool cellsFlag;
	bool bcondFlag;
	bool laplaceFlag;
	bool VariableSizeFlag;
	bool cummulativeVariableSizeFlag;
	bool variableFlag;
	bool refVariableFlag;

// 1. Parameters related to geometrical data (or domain data)
	/// Contains the location of the case
	string path;	
	/// GeometryIO object for reading the data
	Geometry Domain;	
	/// Array of points
	Point* points;	
	/// Array of faces
	Face* faces;	
	/// Array of cells
	Cell* cells;	
	/// Array of boundary conditions
	BoundaryConditions* bcond;	
	/// Number of points, faces, cells, boundaryFaces, internalFaces and boundaryConditions
	int noOfPoints, noOfFaces, noOfCells, noOfBoundaryFaces, noOfBoundaryConditions, noOfInternalFaces;
	/// Reference cells. Refer refCell.cpp and include/cellTypes.h for details 
	RefCell referenceCells[4];
	/// Reference faces.
	RefFace referenceFaces[2];
	/// Shortest distance.
	double shortestDistance;
	/// Characteristic length
	double charLength;
	/// Characteristic speed
	double charSpeed;
	/// CFL number 
	double CFL;

// Constructor and destructor

	/// Constructor function
	DG(string path="/app/Tests/pitzDaily", int order=1, IntFlag::intflag intType = IntFlag::inexact);			

	///destructor function
	~DG(); 

// Other functions
	
	/// Get order of reconstruction 
	int getOrder()const;
	/// Set order of reconstruction
	void setOrder(int order);

	/// Get path of the application
	string getPath()const;
	/// Set path of the application
	void setPath(string path);

	/// Setup reference cells array
	void setupReferenceCellsArray();

	/// Get pointer to reference cells array
	RefCell* getReferenceCellsArray();

	/// Setup reference faces array
	void setupReferenceFacesArray();

	/// Get pointer to reference faces array
	RefFace* getReferenceFacesArray();

	/// Set number of variables
	void setNoOfVariable(int noOfVariable);

	/// Get number of variables
	int getNoOfVariable();

	/// Setter function for name of the dependent variables
	void setVariableName(string VariableName[]);	

	/// Setter function for the size of the variables
	void setVariableSize(int VariableSize[]);	

	/// Set functional details
	void setFunctionalDetails(int order, IntFlag::intflag intType, Solver::solver sol, System::system sys); 

	/// Assign system of equation. Sets the internal flux function.
	void assignSystemOfEquations();

	/// A function pointer for finding the internal flux. Points to the appropriate flux function specified in InternalFlux.h depending on the system of equation.
	void (* findInternalFlux)(Cell *cell, int rkstep); 

	/// A function pointer for solving the Riemann solver. Points to the appropriate Riemann solver specified in RiemannSolvers.h depending on the system of equation and Riemann solver choice. 
	void (* solveRiemannProblem)(Face *face, int rkstep); 

	/// Assigns the appropriate Riemann solver.
	void assignRiemannSolver();

	/// Print all details
	void printFunctionalDetails();			

	/// Set the time step for the simulation
	void setDeltaT(double deltaT);			
        /// Get the time step for the simulation
	double getDeltaT();                              

	/// Set the starting time
	void setStartTime(double startTime);		
        /// Get the starting time
	double getStartTime();                         

	/// Set the ending time
	void setEndTime(double endTime);		
        /// Get the ending time
	double getEndTime();                            
	
	/// Set the time step after which data dumping is done
	void setPrintTime(double printTime);		
        /// Get the time step after which data dumping is done
	double getPrintTime();                          

	/// Set the Integrator Type
	void setIntegratorType(int integratorType);	
        /// Get the Integrator Type in string format
	string getIntegratorType();                     

	/// Get the numerical quadrature type (exact or inexact integration)
	string getQuadratureType();			

	/// Total number of variables requiring time integration
	void setTotNoOfVariableTimeInt(int totNoOfVariableTimeInt);

	/// Get total number of variable
	int getTotNoOfVariables();

	/// Input all details regarding time step, time and integrator in a single function
	void setTemporalDetails(double deltaT, double startTime, double endTime, double printTime, int integratorType);						
	/// Print all temporal details
	void printTemporalDetails();			

	/// For calculating the total number of dependent variables if all were scalar quantitites
	void addTotNoOfVariable();			

	// Allocating space to the 2D variable matrix
	void allocateVariableArraySize();		

	/// Allocating space to the 2D variable matrix in boundary Element
	void allocateBoundaryVariableArraySize();	

	/// File read from the given folder in OpenFOAM format. Default 0.
	void readVariableArray(double Time = 0);	
	
	/// File written to Time folder in OpenFOAM format
	void writeVariableArray(double Time);		

	/// File written to Time folder in OpenFOAM format
	void writeVariableRefArray(double Time);		

	/// Assign variable to cell DOF locations
	void assignVariabletoCell();		

	/// Integrate the cell variable over all DOF locations to get the average value of the cell variable
	void integrateCellVariable();		

	/// A flux vector is configured for each face. 
	void assignFluxVectortoFace();		

	/// Create domain 
	/// Initializes arrays and fills data from the case file
	void createDomain(); 				

	/// Prints the Geometry features 
	void printDomain();	 			

	/// Saves global (x,y,z) location for all the quadrature points for all the cells (3D)
	void getQuadPointsGlobalLocation();

	/// Saves global (x,y,z) location for all the DOF points for all the cells (3D)
	void getDOFPointsGlobalLocation();

	/// Calculate det(J) array for each cell. Stored in Cell::J 
	void calculateJacobian();

	/// Calculate the inverse Jacobian (metric of transformation or dRST_by_dXYZ)
	void calculateInverseJacobian();

	/// Get global (x,y,z) location for all the quadrature points for all the faces (2D)
	void getFaceQuadPointsGlobalLocation(); 

	/// Saves global (x,y,z) location for all the DOF points for all the faces (2D)
	void getFaceDOFPointsGlobalLocation(); 

	/// Maps cell DOF points to face quadrature points (fills array mapOwner/NeighbourDOFPoints)
	void mapFaceDOFPointsToCellDOFPoints();

	/// Maps cell quadrature points to face quadrature points (fills array mapOwner/NeighbourQuadPoints)
	void mapFaceQuadPointsToCellQuadPoints();

	/// Calculate det(J) array for all the quadrature points on the face
	void calculateFaceJacobian(); 

	/// Get shortest distance 
	double getShortestDistance();

	/// Get the characteristic length
	double getCharacteristicLength();

	/// Assign boundary condition functions (refer BoundaryConditions::func2bcond function pointer)
	void assignBoundaryConditions();

	/// Takes variable and size as input and generates space for variable vector
	void generateVariableArray(int noOfVariable, string VariableName[], int VariableSize[]);

	/// Assign a name for each face. If the face is one of the internal faces, then the name is "internal" else it is one of the boundary names
	void assignFaceName();

 	/// Sets the 'noOfBoundaryFaces' in each cell. 0 for internal cells.
	void assignNoOfBoundaryFacesToCells();

	/// Reads variable arrays from folder 'Time'
	void readData(double Time = 0.0);

 	/// Apply the initial conditions
	void applyIC();

	// Depends on the problem. For problems with IC as the reference solution, nothing is done 
	// since the variableRef array already stores IC. 
	// For problems with analytical solution avaialble, this is set to compute the error norms
	void computeReferenceSolution(double time); 

	/// Computes L_infty, L_1 and L_2 errors between the solution and the variableRef array
	void computeError(); 	

	/// Get error array
	double* getError();

	/// Array of errors for all the cells
	double error[3]; 				

	/// Error history over time
	TensorO2<double> ERROR;

	/// Computes M,D,F,S matrices (and anything else that may be needed) for all cells
	void computeCellMatrix();

	/// Computes internal and Riemann fluxes at all the DOF locations
	void computeFlux(int rkstep);

	/// Computes the boundary conditions for all the faces
	void computeBoundaryConditions(int rkstep);
	
	/// Compute the flux residual (Res)
	void computeRES(int rkStep);

	/// Run the application (the binding function)
	void runApplication();

	/// Calculate the characteristic speed (for computation of CFL number)
	void calculateCFL();

	/// RK3 timestepping
	void integratorRK3();



private:
// Variables
// 1. Control parameters
	/// Order of the Interpolation method used
	int order;						
	/// Integration type: 0=inexact, 1=exact integration.
	IntFlag::intflag intType; 						
	/// Which Riemann solver. Refer include/solvers.h
	Solver::solver solverType;  					
	/// Which system of equations. Refer include/systems.h
	System::system systemType;				

// 2. Parameters related to conserved (or primitive) variables
	/// No of dependent variables
	int noOfVariable, totNoOfVariable, totNoOfVariableTimeInt;				
	/// Name of dependent variables: Also serves as input files in folder 0
	string* VariableName; 					
	/// Variablesize = 1 for scalar, Variablesize > 1 for vector
	TensorO1<int> VariableSize;
	/// Cummulative variable size (including vector valued variables)
	TensorO1<int> cummulativeVariableSize;
	/// Saves all data of all dependent variables at cell centres
	TensorO2<double> variable;
	/// reference variable =   analytical value (to compute numerical order)
	//                        or previous time-step value (to compute residue)
	TensorO2<double> variableRef;				
	// Analytical Sod's solution
	TensorO2<double> SodAnalytical;


// 3. Temporal parameters
	/// Starting time
	double startTime;				
	/// Ending time
	double endTime;				
	/// Delta time, timestep
	double deltaT;					
	/// Time after which printing data is to be done
	double printTime;				
	/// Integrator 0: Euler, 1-4: Runge Kutta
	int integratorType;										

};


#endif
