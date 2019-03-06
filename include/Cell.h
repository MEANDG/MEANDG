#ifndef CELL_H
#define CELL_H

#include "sysInclude.h"
#include "Point.h"
#include "Face.h"
#include "mathFunctions.h"
#include "RefCell.h"

class RefCell;

class Face;

class Cell {
public:
	/// Constructor
	Cell (int Id=0);

	/// Destructor
	~Cell();

	bool JacobianFlag;
	bool refCellFlag;

	/// Set unique ID of the cell 
	void setId(int Id);

	/// Get unique ID of the cell 
	int getId();

	/// Get total number of defining points of the cell
	unsigned int getNoOfPoints() const;

	/// Get a particular (ith) defining point of the cell
	Point* getDefiningPoint(int i) const;

	/// Add a defining point to the array of defining points
	void addDefiningPoint(Point *definingPoint);

	/// Set the cell type based on the total number of points
	void setCellType();

	/// Get the cell type based on the total number of points
	CellType::cellType getCellType()const;

	/// Set number of defining faces for the cell
	void setNoOfFaces(unsigned int noOfFaces);

	/// Get number of defining faces for the cell
	unsigned int getNoOfFaces()const;

	/// Add a defining face to the array of defining faces
	void addDefiningFace(Face *definingFace);

	/// Get a defining face from the array of defining faces
	Face * getDefiningFace(int i);

	/// This function checks if the argument point already exists in the array of definingPoints
	bool findPoint(Point *p);

	/// This function finds all the defining points from the array of definingFaces. 
	void findDefiningPoints();

	/// This function returns ordered point index
	int getOrderedPoints(int i)const;

	/// This function returns ordered point index
	int getOrderedFaces(int i);

	/// This array stores the relationship of the current cell with each one of its faces
	/// E.g. faceRelationshipArray = [0, 0, 0, 0, 1, 1] means, the cell is owner (0) for first four faces and neighbour for the last two.
	TensorO1<int> faceRelationshipArray;

	/// Set the largest face
	void setLargestFace(Face *face);

	/// Get the largest face
	Face* getLargestFace();

	/// Get the number of boundary faces
	int getNoOfBoundaryFaces();

	/// Set the number of the boundary faces
	void setNoOfBoundaryFaces(int NBF);
	
	/// Adds the NBF number of faces to existing noOfBoundaryFaces for this cell
	void addNoOfBoundaryFaces(int NBF);

	/// Set the center of the cell
	void calculateCenter();

	/// Get the center of the cell
	Point* getCenter();

	/// Get the pointer to the reference cell corresponding to this cell.
	RefCell* getReferenceCell();

	/// Set the pointer to the reference cell corresponding to this cell.
	void setReferenceCell(RefCell* refCell);

	/// Get degrees of freedom number
	unsigned int getNDOF()const;

	/// Set DOF number
	void setNDOF(unsigned int nDOF);

	/// Set number of quadrature points in the cell
	void setNQuad(unsigned int nQuad);
	
	/// Get number of quadrature points in the cell
	unsigned int getNQuad()const;

	/// Get the pointer to the Jacobian vector (det(J))
	TensorO1<double>* getJacobian();

	/// Get the pointer to the inverse Jacobian vector (at all quad points)
	TensorO3<double>* getInverseJacobian();

	/// Get global location of quadrature points (pointer to TensorO2<> array)
	TensorO2<double>* getQuadPointsGlobalLocation();

	/// Get global location of DOF points (pointer to TensorO2<> array)
	TensorO2<double>* getDOFPointsGlobalLocation();

	/// Map the given local (r,s,t) location to global (x,y,z) coordinate system.
	void mapRSTtoXYZ3DTensor(Point* localPt, Point* globalPt);

	/// Get vertex sign  for particular vertex and for particular dimension
	double getVertexSign(int vertexNo, int coord);

	/// Get entire vertexSign array
	TensorO2<double>* getVertexSignArray();

	/// Assert that the vertex signs are correctly assigned
	void assertVertexSigns();

	/// Set vertex sign 
	void setVertexSign(int vertexNo, int coord, double sign);

	/// Get det(J) value for any given local coordinate point (r,s,t)
	double calculateJacobian3DTensor(double r, double s, double t);

	/// Calculate Inverse Jacobian matrix for given (r,s,t) quadrature point.
	double calculateInverseJacobianMatrix3DTensor(double r,double s,double t, Matrix<double> *dRST_by_dXYZ);

	/// Get pointer to the Mass Matrix
	Matrix<double>* getMassMatrix();

	/// Get pointer to the inverse of the Mass Matrix
	Matrix<double>* getInverseMassMatrix();

	/// Get Vandermonde matrix
	Matrix<double>* getVandermondeMatrix();

	/// Get Dx matrix
	Matrix<double>* getDxMatrix();

	/// Get Dy matrix
	Matrix<double>* getDyMatrix();

	/// Get Dz matrix
	Matrix<double>* getDzMatrix();

	/// Get Dtildex matrix
	Matrix<double>* getDtildexMatrix();

	/// Get Dtildey matrix
	Matrix<double>* getDtildeyMatrix();

	/// Get Dtildez matrix
	Matrix<double>* getDtildezMatrix();

	/// Get flux array
	TensorO3<double>* getFMatrix();

	/// Get flux array, multiplied by MI
	TensorO3<double>* getFtildeMatrix();

	/// Get the filter matrix
	Matrix<double>* getFilterMatrix();

	/// Resize the matrices
	void resizeMatrices(unsigned int Size);

	/// Get Variable vector
	TensorO3<double>* getVariableArray();

	/// Get the reference variable array
	TensorO2<double>* getVariableRefArray();

	/// Get a particular variable value
	double getVariable(int rkStep, int varNo, int DOF);

	/// Get Variable residual
	TensorO3<double>* getVariableResidualArray();

	/// Get analytical flux vector in x direction (i.e. fluxVectorx)
	TensorO2<double>* getFluxVectorx();

	/// Get analytical flux vector in y direction (i.e. fluxVectory)
	TensorO2<double>* getFluxVectory();

	/// Get analytical flux vector in z direction (i.e. fluxVectorz)
	TensorO2<double>* getFluxVectorz();

	/// Get flux star (Riemann) flux vector for the cell (for all the faces)
	TensorO3<double>* getFluxStarVector();

	// Printing Functions
	/// Print Mass Matrix on screen
	void printMassMatrix();
	/// Print Mass Matrix in a file
	void printMassMatrix(ofstream &file);

	/// Print Inverse of the mass matrix
	void printInverseMassMatrix();
	void printInverseMassMatrix(ofstream &file);

	/// Print Vandermonde matrix
	void printVandermondeMatrix();
	void printVandermondeMatrix(ofstream &file);

	/// Print differentiation matrix (Dx, Dy, Dz)
	void printDifferentialMatrix();
	void printDifferentialMatrix(ofstream &file);

	/// Print differentiation matrix multiplied with M_inv (Dtildex, Dtildey, Dtildez)
	void printDTildeMatrix();
	void printDTildeMatrix(ofstream &file);

	/// Print the flux matrix 
	void printFluxMatrix();
	void printFluxMatrix(ofstream &file);

	/// Print the Ftilde matrix
	void printFTildeMatrix();
	void printFTildeMatrix(ofstream &file);

	/// Print the filter matrix
	void printFilterMatrix();
	void printFilterMatrix(ofstream &file);

private:

	/// A unique id for each cell. Set during initialization.
	unsigned int id;

	/// Number of defining points 
	unsigned int noOfPoints;

	/// Array of defining points
	Point* definingPoints[8];

	/// Ordered sequence of the defining points
	int orderedPoints[8];

	/// Ordered sequence of the defining faces
	int orderedFaces[6];

	/// Cell type (Hex, Tet, Prism or Pyramid). Refer include/cellTypes.h
	CellType::cellType celltype;
	bool cellTypeFlag;		//no further point can be added once the flag is set to true. Default false.

	/// Number of defining faces
	unsigned int noOfFaces;

	/// Array of defining faces
	Face* definingFaces[8];

	/// Largest face of the cell
	Face* largestFace;

	/// Number of boundary faces that make up the given cell, 0 if none of the faces is at the boundary of the domain
	int noOfBoundaryFaces;

	// Function to order the defining points in the correct sequence
	void orderDefiningPoints();

	/// Center point of the cell
	Point center;

	/// Reference cell corresponding to this cell
	RefCell* refCell;

	/// Number of degrees of freedom
	unsigned int nDOF;
	bool nDOFFlag;

	/// Number of quadrature points
	unsigned int nQuad;
	bool nQuadFlag;

	/// Jacobian vector (det(J) actually)
	TensorO1<double> J;

	/// Vertex signs (for shape function calculation). Dim: noOfPoints x 3  (for x, y, z)
	TensorO2<double> vertexSign;

	/// Inverse Jacobian matrix or metric terms
	TensorO3<double> dRST_by_dXYZ;			// dim: nQuad x 3 x 3 i.e. Inverse Jacobian at nQuad points
	
	/// Quadrature points global location
	TensorO2<double> quadPointsGlobalLocation; 	//dim: nQuad x 3, stores x,y,z location of all the quadrature points.
                                                        //implementation is done in DG::getQuadPointsLocation();

	/// DOF points global location 
	TensorO2<double> DOFPointsGlobalLocation; 	//dim: nDOF x 3, stores x,y,z location of all the DOF points.
							//implementation is done in DG::getDOFPointsLocation();

	/// Mass Matrix
	Matrix<double> M;

	/// Inverse of the Mass Matrix
	Matrix<double> MI;

	/// Vandermonde Matrix
	Matrix<double> V;

	/// Dx (differentiation or stiffness (depending on the notation)) matrix with gradient in x direction
	Matrix<double> Dx;

	/// Dy (differentiation or stiffness (depending on the notation)) matrix with gradient in y direction
	Matrix<double> Dy;

	/// Dz (differentiation or stiffness (depending on the notation)) matrix with gradient in z direction
	Matrix<double> Dz;

	/// Dtildex (differentiation or stiffness (depending on the notation)) matrix with gradient in x direction, multiplied with MI
	Matrix<double> Dtildex;

	/// Dtildey (differentiation or stiffness (depending on the notation)) matrix with gradient in y direction, multiplied with MI
	Matrix<double> Dtildey;

	/// Dtildez (differentiation or stiffness (depending on the notation)) matrix with gradient in z direction, multiplied with MI
	Matrix<double> Dtildez;

	/// F (flux) matrix
	TensorO3<double> F;

	/// Ftilde (flux) matrix, multiplied by MI
	TensorO3<double> Ftilde;

	/// Filter Matrix
	Matrix<double> Filter;


	/// Matrix that stores the values of properties at all DOF points. 
	/// First index: RK step,   Second index: property type (e.g. density),     Third index: DOF point number
	TensorO3<double> variable;

	/// Reference variable at the DOF location
	/// First index: property type(e.g. density),  Second index: DOF point number
	TensorO2<double> variableRef;

	/// Matrix that stores the values of variable residual at all DOF points. 
	/// First index: RK step,   Second index: perperty type (e.g. density),     Third index: DOF point number
	TensorO3<double> variableResidual;

	/// Matrix that stores the analytical value of the flux for all the DOF locations.
	/// First index: variable number (e.g. density),   Second index: DOF location number
	/// To be multiplied by the differentiation matrix later
	TensorO2<double> fluxVectorx;
	TensorO2<double> fluxVectory;
	TensorO2<double> fluxVectorz;

	/// The array storing the Riemann flux found out at all the faces. These are assembled in one global vector for the cell.
	/// This is to be mulitplied by the flux matrix later.
	/// First index: Face number,    Second index: RK step,    Third index: DOF location number
	TensorO3<double> fluxStarVector;

};


namespace Test{
	void TestCell();
	void TestCellInitialization();
	void TestDefiningPointsForCell();
	void TestDefiningFacesForCell();
};

#endif
