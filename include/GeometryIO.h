#ifndef GEOMETRYIO_H
#define GEOMETRYIO_H

#include "sysInclude.h"
#include "mathFunctions.h"
#include "Point.h"
#include "Face.h"
#include "Cell.h"
#include "BoundaryConditions.h"


class Geometry{
	public:

		/// Constructor
		Geometry();

		/// Destructor
		~Geometry();

		/// Get number of points (as well as implicitly store the value)
		int setNoOfPoints(string &dir);

		/// Get number of faces 
		int setNoOfFaces(string &dir);

		/// Get number of cells
		int setNoOfCells(string &dir);

		/// Get number of boundary conditions
		int setNoOfBoundaryConditions(string &dir);

		/// Get number of boundary faces
		int setNoOfBoundaryFaces(string &dir);

		/// Run all the getFunctions to get the size of the problem
		void getProblemData(string &dir);


		/// Get number of points (as well as implicitly store the value)
		int getNoOfPoints()const;

		/// Get number of faces 
		int getNoOfFaces()const;

		/// Get number of cells
		int getNoOfCells()const;

		/// Get number of boundary conditions
		int getNoOfBoundaryConditions()const;

		/// Get number of boundary faces
		int getNoOfBoundaryFaces()const;

		/// Get number of boundary faces
		int getNoOfInternalFaces()const;

		/// Function to read the points file in OpenFOAM format 
		/// takes the name of point file and a pointer to array of objects of Point class as an input
		/// The function reads the coordinates of the points and adds them to the list vector
		void readPointFile(string& pointFile, Point *p);

		/// Function to read the face file in OpenFOAM format 
		/// takes the name of face file and pointers to arrays of objects of Face class and Point class as an input
		/// The function reads the points that contsruct a face and add the face to the list vector
		void readFaceFile(string& faceFile, Face* faces, Point* points);

		/// Function to read the boundary conditions files in OpenFOAM format
		/// takes the name of the boundary file and appends the array of boundary condition 
		void readBoundaryFile(string& boundaryFile, BoundaryConditions* bcond);

		/// Function to read the owner file in OpenFOAM format 
		/// takes the name of face file and pointers to arrays of objects of Face class and Cell class as an input
		/// The function reads the owners of the faces and simultaneously builds the cells using the faces
		void readOwnerFile(string& ownerFile, Cell* cells, Face* faces);

		/// Function to read the neighbour file in OpenFOAM format 
		/// takes the name of face file and vector of cell class after the owner faces are linked as an input
		/// The function reads the neighbour of the faces and simultaneously builds the cells using the faces
		void appendNeighbourFile(string& neighbourFile, Cell* c, Face* faces);

		/// Function to assign owner and neighbour cells to each of the faces
		/// The faces lying on the boundary of the domain are assigned with a NULL cell. 
		/// Such faces show "GHOST" as the neighbour cell when inquired. 
		void assignOwnerCellsToFaces(string& ownerFile, Cell* cells, Face* faces);
		void assignNeighbourCellsToFaces(string& neighbourFile, Cell* cells, Face* faces);

		/// Following function fills an array in each cell which tells whether a face from definingFaces has the current
		/// cell as owner or neighbour.
		void fillFaceRelationshipArray(Cell* cells, Face* faces);

		/// Sets the face->normalDirFlag (refer include/Face.h) 
		void findFaceNormalsOrientation(Cell* cells, Face* faces);

		/// Find the shortest distance in the domain
		double findShortestDistance(Cell* cells);

		/// Get shortest distance on the domain
		double getShortestDistance();

		/// Find the characteristic length of the domain
		double findCharacteristicLength(Cell* cells);

		/// Get the (precomputed) characteristic length
		double getCharacteristicLength();

		/// Assign vertex signs to all the cells (for Jacobian calculations and shape function calculations) 
		void assignVertexSigns(Cell *cells, Face *faces);

		/// Fill all the data arrays
		void fillDataArrays(string& dir, Point* points, Face* faces, Cell* cells, BoundaryConditions* bcond);

		/// Print the domain details
		void print();

	private:

		/// Number of points in the domain
		int noOfPoints;

		/// Number of faces in the domain
		int noOfFaces;

		/// Number of cells in the domain
		int noOfCells;

		/// Number of boundary faces in the domain
		int noOfBoundaryFaces;

		/// Number of internal faces in the domain
		int noOfInternalFaces;

		/// Number of boundary conditions
		int noOfBoundaryConditions;

		/// Shortest distance on the domain
		double shortestDistance;

		/// Characteristic length
		double charLength;


};

namespace Test{
	void TestGeometryIO();
};

#endif
