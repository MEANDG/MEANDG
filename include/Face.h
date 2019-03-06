#ifndef FACE_H
#define FACE_H

#include "sysInclude.h"
#include "Point.h"
#include "mathFunctions.h"
#include "RefFace.h"

class RefFace;

class Cell; //forward declaration


class Face {

	
public:
	
	/// Constructor with dfault value pointing at the origin
	Face(int Id=0);
	/// Destructor
	~Face();

	/// Mapping flag
	bool mapFlag;
	/// Jacobian flag
	bool JacobianFlag;

	/// Set unique ID of the face 
	void setId(unsigned int Id);

	/// Get unique ID of the face 
	int getId();

	/// Set the name of the face
	void setName(string name);

	/// Get the name of the face
	string getName();

	/// Get total number of defining points of the face
	unsigned int getNoOfPoints() const;

	/// Get a particular (ith) defining point of the face
	Point* getDefiningPoint(int i) const;

	/// Add a defining point to the array of defining points
	void addDefiningPoint(Point *definingPoint);

	/// Set the face type based on the total number of points
	void setFaceType();

	/// Get the face type based on the total number of points
	FaceType::faceType getFaceType()const;

	/// Set the owner cell pointer
	void setOwnerCell(Cell *cell);

	/// Get the owner cell pointer
	Cell* getOwnerCell();

	/// Set the neighbour cell pointer
	void setNeighbourCell(Cell *cell);

	/// Get the neighbour cell pointer
	Cell* getNeighbourCell();

	/// Flag which tells if the given face is a boundary face (default true, set to false if neighbour present);
	bool isBoundary();

	/// Calculate normal vector 
	void calculateNormal();

	/// Get the pointer to the normal vector (of the type TensorO1<double>)
	TensorO1<double>* getNormal();
	TensorO1<double>* getTangent1();
	TensorO1<double>* getTangent2();


	/// Get area of the face
	double getArea()const;

	/// Compute area of the face
	void calculateArea();

	/// Calculate center of the face
	void calculateCenter();
	
	/// Get center of the face
	Point* getCenter();

	/// Get ordered points sequence
	int getOrderedPoints(int index);

	/// Find whether the argument point exists in the face.definingPoints[]
	bool findPoint(Point *p);

	/// Flag indicating whether normal points into the owner cell or away
        /// normalDirFlag = 1.0 -> normal vector points away from the owner cell (as expected), else -1.0
	/// It is a double dtype because it is directly involved in the computations.
	double normalDirFlag;

	/// Get method for numerical flux for nth variable at mth DOF
	double getFlux(int variable, int DOF);

	/// Get flux vector
	TensorO2<double>* getFluxArray();

	/// Set method for numerical flux for nth variable at mth DOF
	void setFlux(int variable, int DOF, double Value);

	/// Get method for numerical flux for nth variable at mth DOF
	double getPrimitiveVariable(int variable, int DOF);

	/// Get pointer to primitiveVariableVector array
	TensorO2<double>* getPrimitiveVariableArray();

	/// Set method for numerical flux for nth variable at mth DOF
	void setPrimitiveVariable(int variable, int DOF, double Value);

	/// Get number of DOFs 
	int getNDOF()const;

	/// Set number of DOFs 
	void setNDOF(int nDOF);

	/// Get number of Quadrature points 
	int getNQuad()const;

	/// Set number of Quadrature points 
	void setNQuad(int nQuad);

	/// Get the pointer to the reference cell corresponding to this cell.
	RefFace* getReferenceFace();

	/// Set the pointer to the reference cell corresponding to this cell.
	void setReferenceFace(RefFace* refFace);

	/// Get global location of quadrature points (pointer to TensorO2<> array)
	TensorO2<double>* getQuadPointsGlobalLocation();

	/// Get global location of DOF points (pointer to TensorO2<> array)
	TensorO2<double>* getDOFPointsGlobalLocation();

	/// Get vertex sign  for particular vertex and for particular dimension
	double getVertexSign(int vertexNo, int coord);

	/// Get entire vertexSign array
	TensorO2<double>* getVertexSignArray();

	/// Assert that the vertex signs are correctly assigned
	void assertVertexSigns();

	/// Set vertex sign 
	void setVertexSign(int vertexNo, int coord, double sign);


	/// Get the pointer to rotationMatrixParallelToXY
	Matrix<double>* getRotationMatrixParallelToXY();

	/// Get the pointer to inverse rotationMatrixParallelToXY
	Matrix<double>* getInverseRotationMatrixParallelToXY();

	/// Rotate the given vector (TensorO1<double>> of size 3) 
	void rotateParallelToXY(TensorO1<double> *OldPt, TensorO1<double> *newPt);

	/// Rotate the given vector (TensorO1<double>> of size 3) back to original direction 
	void rotateBackParallelToXY(TensorO1<double> *OldPt, TensorO1<double> *newPt);

	/// Calculate rotation matrix which rotates the given face || to xy plane (i.e. normal pointing in the +ve z direction)
	void calculateRotationMatrixParallelToXY();

	/// Calculate inverse rotation matrix
	void calculateInverseRotationMatrixParallelToXY();



	/// Get the pointer to rotationMatrixNormal
	Matrix<double>* getRotationMatrixNormal();

	/// Get the pointer to inverse rotationMatrixNormal
	Matrix<double>* getInverseRotationMatrixNormal();

	/// Rotate the given vector (TensorO1<double> of size 3) in normal-tangent1-tangent2 coordinate system 
	void rotateNormal(TensorO1<double> *OldPt, TensorO1<double> *newPt);

	/// Rotate the given vector (TensorO1<double> of size 3) back to original direction 
	void rotateBackNormal(TensorO1<double> *OldPt, TensorO1<double> *newPt);

	/// Calculate rotation matrix projects the given vector in normal-tangent1-tangent2 coordinate system
	void calculateRotationMatrixNormal();

	/// Calculate inverse rotation matrix
	void calculateInverseRotationMatrixNormal();



	/// Map the local coordinate to global coordinate
	void mapRSTtoXYZQuad(Point *localPt, Point *globalPt);

	/// Calculate Jacobian of the quadrilateral face
	double calculateJacobianQuadFace(double r, double s);

	/// Get the pointer to the Jacobian vector (det(J))
	TensorO1<double>* getJacobian();

	/// Get owner quad point corresponding to the arg quadrature point
	int getOwnerQuadPoint(int faceQP);

	/// Get neighbour quad point corresponding to the arg quadrature point
	int getNeighbourQuadPoint(int faceQP);

	/// Get owner DOF point corresponding to the arg DOFrature point
	int getOwnerDOFPoint(int faceDOF);

	/// Get neighbour DOF point corresponding to the arg DOFrature point
	int getNeighbourDOFPoint(int faceDOF);

	/// Get pointer to mapOwnerQuadPoints
	TensorO1<int>* getOwnerQuadPointsArray();

	/// Get pointer to mapNeighbourQuadPoints
	TensorO1<int>* getNeighbourQuadPointsArray();

	/// Get pointer to mapOwnerDOFPoints
	TensorO1<int>* getOwnerDOFPointsArray();

	/// Get pointer to mapNeighbourDOFPoints
	TensorO1<int>* getNeighbourDOFPointsArray();

private:

	/// A unique id for each face. Set during initialization.
	unsigned int id;

	/// A name for the face ("internal" or one of the boundary names)
	string name;

	/// Number of defining points 
	unsigned int noOfPoints;

	/// Array of defining points
	Point* definingPoints[4];

	/// Face center
	Point center;

	/// Face type (tri or quad). Refer include/faceTypes.h
	FaceType::faceType facetype;
	bool faceTypeFlag;		//no further point can be added once the flag is set to true. Default false.


	/// Pointer to the owner cell (object of the class Cell)
	Cell *owner;
	/// Pointer to the neighbour cell (object of the class Cell)
	Cell *neighbour;

	/// Boundary flag. Set to false if neighbour is present. Default true.
	bool boundaryFlag;

	/// Face normal 
	TensorO1<double> normal;
	/// Face tangent 1
	TensorO1<double> tangent1;
	/// Face tangent 2
	TensorO1<double> tangent2;
	
	/// A flag which is turned to true when the normal is computed
	bool normalFlag;


	/// Area of the face
	double area;

	/// Ordered sequence of the defining points
	int orderedPoints[4];


	void orderDefiningPoints();


	/// Numerical flux at the face, rows: conservative variables,  columns: DOF location
	TensorO2<double> Flux;

	/// Primitive variable vector at the face, rows: primitive variables,  columns: DOF location
	TensorO2<double> primitiveVariableVector;


	/// Number of DOF location on the face
	int nDOF;

	// Flag for setting DOF number
	bool nDOFFlag;

	/// Number of the quadrature points on the face
	int nQuad;

	// Flag for setting Quad points number
	bool nQuadFlag;


	/// Vertex signs (for shape function calculation). Dim: noOfPoints x 2  (for x, y)
	TensorO2<double> vertexSign;

	/// Quadrature points global location
	TensorO2<double> quadPointsGlobalLocation; 	//dim: nQuad x 3, stores x,y,z location of all the quadrature points.
                                                        //implementation is done in DG::getQuadPointsLocation();
	/// DOF points global location 
	TensorO2<double> DOFPointsGlobalLocation; 	//dim: nDOF x 3, stores x,y,z location of all the DOF points.
							//implementation is done in DG::getDOFPointsLocation();

	/// Reference face corresponding to this face
	RefFace* refFace;
	bool refFaceFlag;


	/// Matrix which rotates the given face || to xy plane (i.e. normal vector pointing in the +ve z direction).
	Matrix<double> rotationMatrixParallelToXY;

	/// Rotates the face back to its original orientation
	Matrix<double> inverseRotationMatrixParallelToXY; 


	/// Matrix which projects the given vector in normal-tangent1-tangent2 coordinate system
	Matrix<double> rotationMatrixNormal;

	/// Rotates the face back to its original orientation
	Matrix<double> inverseRotationMatrixNormal; 


	/// Jacobian vector (det(J) actually)
	TensorO1<double> J;

	/// This array maps the 2D quadrature point for surface integration to 3D global quadrature point
        /// i.e. which 3D global quadrature point is same as the local face quadrature point.
        /// Just gives a single array of global quad point numbers. It is understood that their index corresponds to the local number.
	/// e.g. array [3,4,2,1] means the face quadrature point 0 = global quadrature point 3, face q p 1 = global q p 4 etc.
	/// This array will be used for constructing the LagMatrix for the surface integration
	TensorO1<int> mapOwnerQuadPoints;
	TensorO1<int> mapNeighbourQuadPoints;

	/// Similar to mapOwnerQuadPoints, but for DOF locations
	/// For inexact integration, both mapOwnerQuadPoints and mapOwnerDOFPoints arrays are same (same for neighbour)
	TensorO1<int> mapOwnerDOFPoints;
	TensorO1<int> mapNeighbourDOFPoints;




	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Prevent copy constructor
	Face(const Face&);
	// copy assignment operator 
	Face& operator=(const Face& tmp_obj){ 
	        return *this; 
	};


};


/********************************Tests***********************************/

namespace Test{
	void TestFace();
	void TestFaceInitialization();
	void TestDefiningPoints();
	void TestNormalAndTangents();
	void TestOwnerAndNeighbourCells();
	void TestFaceArea();
	void TestCenter();
};
#endif
