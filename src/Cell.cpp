#include "Cell.h"


Cell::Cell(int Id){
	this->id = Id;
	this->noOfPoints = 0;
	this->noOfFaces = 0;
	// Default order of the defining points
	for (Index i=0; i<8; i++){
		this->orderedPoints[i] = i;
	};

	// Default order of the defining faces
	for (Index i=0; i<6; i++){
		this->orderedFaces[i] = i;
	};

	this->cellTypeFlag = false;
	this->faceRelationshipArray.setSize(6);

	this->nDOFFlag = false;
	this->nQuadFlag = false;
	this->refCellFlag = false;
	this->JacobianFlag = false;

};


Cell::~Cell(){
};


void Cell::setId(int Id){
	this->id = Id;
};

int Cell::getId(){
	return this->id;
};

unsigned int Cell::getNoOfPoints() const{
	assert(this->noOfPoints >= 4 && "Cell::getNoOfPoints returns error. Make sure that all the points are added\n.");
	return this->noOfPoints;
};

Point* Cell::getDefiningPoint(int i) const{
	assert(i < this->noOfPoints && "Cell::getDefiningPoint() returns error. Argument index greater than total number of defining points.\n");
	return this->definingPoints[this->getOrderedPoints(i)];
};

int Cell::getOrderedPoints(int i)const{
	return this->orderedPoints[i];
};

void Cell::addDefiningPoint(Point *definingPoint){
	assert(this->noOfPoints <= 8 && "Cell::addDefiningPoint() returns error. Cell already has atleast 8 points.\n");
	assert(this->cellTypeFlag == false && "Cell::addDefiningPoint() returns error. The cell type is already finalized.\n");
	this->definingPoints[this->noOfPoints] = definingPoint;
	this->noOfPoints ++;
};

void Cell::setCellType(){
	assert(this->cellTypeFlag == false && "Cell::setCellType() returns error. The cell type is already finalized.\n");
	assert(this->noOfPoints == 4 or this->noOfPoints == 5 or this->noOfPoints == 6  or this->noOfPoints == 8  && "Cell::setCellType() returns error. Build the cell completely first.\n");

	if (this->noOfPoints == 4){
		celltype = CellType::Tet;
	}
	else if (this->noOfPoints == 5){
		celltype = CellType::Pyramid;
	}
	else if (this->noOfPoints == 6){
		celltype = CellType::Prism;
	}
	else if (this->noOfPoints == 8){
		celltype = CellType::Hex;
	};
	this->cellTypeFlag = true; 	//once this flag is set to true, cell type can not be reset or new point can not be added.

	this->orderDefiningPoints();
};


CellType::cellType Cell::getCellType() const{
	assert(this->cellTypeFlag == true && "Cell::setCellType() returns error. Set the celltype first.\n");
	return this->celltype;
};



void Cell::setNoOfFaces(unsigned int noOfFaces){
	assert(noOfFaces == 4 or noOfFaces == 5 or noOfFaces == 6 or noOfFaces == 8 && "Cell::setNoOfFaces() returns error. Please make sure that the argument lies in [4,5,6,8].\n");
	this->noOfFaces = noOfFaces;
};

unsigned int Cell::getNoOfFaces() const{
	assert(noOfFaces > 0 && "Cell::getNoOfFaces() returns error. Make sure that all the faces have been added.\n");
	return this->noOfFaces;
};


void Cell::addDefiningFace(Face *definingFace){
	assert(this->noOfFaces <= 8 && "Cell::addDefiningFace() returns error. Cell already has atleast 8 faces.\n");
	assert(this->cellTypeFlag == false && "Cell::addDefiningFace() returns error. The cell type is already finalized.\n");
	this->definingFaces[this->noOfFaces] = definingFace;
	this->noOfFaces ++;
};

Face* Cell::getDefiningFace(int i){
	assert(i < this->noOfFaces && "Cell::getDefiningFace() returns error. Argument index greater than total number of defining faces.\n");
	return this->definingFaces[this->getOrderedFaces(i)];
};

int Cell::getOrderedFaces(int i){
	return this->orderedFaces[i];
};

bool Cell::findPoint(Point *p){
	bool flag = false;
	assert(this->noOfPoints > 0 && "Cell::findPoint() returns error. The cell has no points right now.\n");
	for (Index i=0; i<this->noOfPoints; i++){
		if (p == this->definingPoints[i]){
			flag = true;
			break;
		}
		else{
			flag = false;
		};
	};

	// If flag == true, then the point already exists else not
	return flag;
};


/// This function fills the array of defining points from the array of pre-filled definingFaces.
void Cell::findDefiningPoints(){
	assert(this->noOfFaces > 3 && "Cell::findDefiningFaces() returns error. Min 4 faces required for a functional cell.\n");

	// First add all the points of the zeroth face to the cell defining Points (which as of now is an empty array)
	for (Index i=0; i < this->definingFaces[0]->getNoOfPoints() ; i++){
		this->addDefiningPoint(this->definingFaces[0]->getDefiningPoint(i));
	};

	// From the face #1 onwards, add only those points which are not already added
	for (Index i=1; i< this->getNoOfFaces() ; i++){
		for (int j=0; j < this->definingFaces[i]->getNoOfPoints() ; j++){
			// the following function checks whether the given point is already added or not
			bool a = this->findPoint(this->definingFaces[i]->getDefiningPoint(j));
			if (a == false){
				// this means point is not already added. Then add it.
				this->addDefiningPoint(this->definingFaces[i]->getDefiningPoint(j));
			};
		};
	};
};


void Cell::orderDefiningPoints(){
	assert(this->cellTypeFlag == true && "Cell::orderDefiningPoints() returns error. Make sure that all the points have been added.\n");

	// In the current code, first the faces are added (based on the OpenFOAM files)
	// Therefore first the ordered array of the faces is obtained.
	// Then Cell::findDefiningPoints() is called to append array of defining points. 

	assert(this->getCellType() == CellType::Hex && "Cell::orderDefiningPoints() returns error. The code is developed only for hex cells as of now.\n");

	if (this->getCellType() == CellType::Hex){

		//_______________________________________________________________________________________________________________
		// First, we find the order of the defining faces.
		// The sequence of the faces follows the same logic as the extruded faces in Gmsh
		// i.e. the base face is 0, followed by the opposite face and then the rest of the faces in the same sequence as the defining points of the base face. 
			// The final sequence of the faces should be:
			// Base: 0,  (containing points 0 1 2 3)
			// Top : 1,  (containing points 4 5 6 7), with 4 opposite of 0, 5 opposite of 1 etc.
			// Front:2,  (sharing pts 0 and 1 with base) 
			// LHS:  3,  (sharing pts 1 and 2 with base)
			// Back: 4,  (sharing pts 2 and 3 with base)
			// RHS:  5,  (sharing pts 3 and 4 with base)
		/*

		    	// Figure: Face orientation and numbering

			                          
				_______________________
				*		       *                             	
				|  *            1      |  *
				|     *      (top)     |     *     
				|        * ____________|________*
				|         |            |          |
				|         |            |          |
				|   2     |  3 (LHS)   |    4     |
		    ----->	| (front) |            |  (back)  |
		(observer	|         |            |          |
			     	|_________|__________5(RHS)       |
				*	  |      	*         |                    	
				   *      |    0(base)     *      |
				      *   |                   *   |
					 *|______________________*|
					  
					  			    

		*/


		// NOTE: this->orderedFaces[0] = 0 by default
		// Find the face opposite to 0th face.

		bool flag = false;
		int faceNo = 0;

		for (Index i=1; i<this->getNoOfFaces(); i++){

			for (Index j=0; j<this->getDefiningFace(0)->getNoOfPoints(); j++){

				Point *pt = this->getDefiningFace(0)->getDefiningPoint(j);

				flag = this->definingFaces[i]->findPoint(pt);

				if (flag == true){
					break;
				};
			};

			if (flag == false){
				faceNo = i;
			};
		};

		this->orderedFaces[1] = faceNo;

		//_____________________________________________________________________
		// Find the face for which first two points of the 0th face are shared

		bool flag1 = false;
		bool flag2 = false;

		Point *pt1 =  this->getDefiningFace(0)->getDefiningPoint(0);
		Point *pt2 =  this->getDefiningFace(0)->getDefiningPoint(1);

		faceNo = 0;
		for (Index i=1; i<this->getNoOfFaces(); i++){

			flag1 = this->definingFaces[i]->findPoint(pt1);
			flag2 = this->definingFaces[i]->findPoint(pt2);

			if (flag1 == true and flag2 == true){
				faceNo = i;
				break;
			};
		};

		this->orderedFaces[2] = faceNo;

		//_____________________________________________________________________
		// Find the rest of the faces using similar procedure. Face: 3

		flag1 = false;
		flag2 = false;

		pt1 =  this->getDefiningFace(0)->getDefiningPoint(1);
		pt2 =  this->getDefiningFace(0)->getDefiningPoint(2);

		faceNo = 0;
		for (Index i=1; i<this->getNoOfFaces(); i++){

			flag1 = this->definingFaces[i]->findPoint(pt1);
			flag2 = this->definingFaces[i]->findPoint(pt2);

			if (flag1 == true and flag2 == true){
				faceNo = i;
				break;
			};
		};

		this->orderedFaces[3] = faceNo;

		//_____________________________________________________________________
		// Find the rest of the faces using similar procedure. Face: 4

		flag1 = false;
		flag2 = false;

		pt1 =  this->getDefiningFace(0)->getDefiningPoint(2);
		pt2 =  this->getDefiningFace(0)->getDefiningPoint(3);

		faceNo = 0;
		for (Index i=1; i<this->getNoOfFaces(); i++){

			flag1 = this->definingFaces[i]->findPoint(pt1);
			flag2 = this->definingFaces[i]->findPoint(pt2);

			if (flag1 == true and flag2 == true){
				faceNo = i;
				break;
			};
		};

		this->orderedFaces[4] = faceNo;

		//_____________________________________________________________________
		// Find the rest of the faces using similar procedure. Face: 5

		flag1 = false;
		flag2 = false;

		pt1 =  this->getDefiningFace(0)->getDefiningPoint(3);
		pt2 =  this->getDefiningFace(0)->getDefiningPoint(0);

		faceNo = 0;
		for (Index i=1; i<this->getNoOfFaces(); i++){

			flag1 = this->definingFaces[i]->findPoint(pt1);
			flag2 = this->definingFaces[i]->findPoint(pt2);

			if (flag1 == true and flag2 == true){
				faceNo = i;
				break;
			};
		};

		this->orderedFaces[5] = faceNo;
		

		//_____________________________________________________________________
		// Now the faces have been ordered properly. Next, we order the defining points.

		/*

		    	// Figure: Vertices orientation and numbering

			
			      5                         6
				_______________________
				*		       *                             	
				|  *                   |  *
				|     *    4           |     *     7
				|        * ____________|________*
				|         |            |          |
				|         |            |          |
				|         |            |          |
				|         |            |          |
				|         |            |          |
			     1	|_________|____________|  2       |
				*	  |      	*         |                    	
				   *      |                *      |
				      *   |                   *   |
					 *|______________________*|
					  
					  0			    3

		*/

		//_______________________________________________________________________________________________________________


		// The first four points are in a sequence corresponding to the 0th face.

		for(Index i=0; i<this->getDefiningFace(0)->getNoOfPoints(); i++){
			this->orderedPoints[i] = i;
		};

		// Next four points are the ones which are directly opposite to the first four on the opposite face of the 0th face


		//_____________________________________________________________________
		// Point 4 is common between faces (ordered) as 1, 2 and 5 (i.e. opposite and two neighbours sharing point 0)
		// Remember, the orderedFaces[] array is finalized. So it will return the faces in the correct order after this point.
		
		for (Index i=0; i<this->getNoOfPoints(); i++){
			Point *pt = this->definingPoints[i];
			
			bool flag1 = false;
			bool flag2 = false;
			bool flag3 = false;

			flag1 = this->getDefiningFace(1)->findPoint(pt);
			flag2 = this->getDefiningFace(2)->findPoint(pt);
			flag3 = this->getDefiningFace(5)->findPoint(pt);

			if (flag1 == true and flag2 == true and flag3 == true){
				this->orderedPoints[4] = i;
			};
		};
		//_____________________________________________________________________
		// Point 5 is common between faces (ordered) as 1, 2 and 3 (i.e. opposite and two neighbours sharing point 1)
		// Remember, the orderedFaces[] array is finalized. So it will return the faces in the correct order after this point.
		
		for (Index i=0; i<this->getNoOfPoints(); i++){
			Point *pt = this->definingPoints[i];
			
			bool flag1 = false;
			bool flag2 = false;
			bool flag3 = false;

			flag1 = this->getDefiningFace(1)->findPoint(pt);
			flag2 = this->getDefiningFace(2)->findPoint(pt);
			flag3 = this->getDefiningFace(3)->findPoint(pt);

			if (flag1 == true and flag2 == true and flag3 == true){
				this->orderedPoints[5] = i;
			};
		};
		//_____________________________________________________________________
		// Point 6 is common between faces (ordered) as 1, 3 and 4 (i.e. opposite and two neighbours sharing point 2)
		// Remember, the orderedFaces[] array is finalized. So it will return the faces in the correct order after this point.
		
		for (Index i=0; i<this->getNoOfPoints(); i++){
			Point *pt = this->definingPoints[i];
			
			bool flag1 = false;
			bool flag2 = false;
			bool flag3 = false;

			flag1 = this->getDefiningFace(1)->findPoint(pt);
			flag2 = this->getDefiningFace(3)->findPoint(pt);
			flag3 = this->getDefiningFace(4)->findPoint(pt);

			if (flag1 == true and flag2 == true and flag3 == true){
				this->orderedPoints[6] = i;
			};
		};
		//_____________________________________________________________________
		// Point 7 is common between faces (ordered) as 1, 4 and 5 (i.e. opposite and two neighbours sharing point 3)
		// Remember, the orderedFaces[] array is finalized. So it will return the faces in the correct order after this point.
		
		for (Index i=0; i<this->getNoOfPoints(); i++){
			Point *pt = this->definingPoints[i];
			
			bool flag1 = false;
			bool flag2 = false;
			bool flag3 = false;

			flag1 = this->getDefiningFace(1)->findPoint(pt);
			flag2 = this->getDefiningFace(4)->findPoint(pt);
			flag3 = this->getDefiningFace(5)->findPoint(pt);

			if (flag1 == true and flag2 == true and flag3 == true){
				this->orderedPoints[7] = i;
			};
		};

		// Also, we need to assign appropriate space for Cell::vertexSign array. 
		// This array will be used in the shape function calculation as well as Jacobian calculation
		this->vertexSign.setSize(this->getNoOfPoints(), 3);	// first index on defining points, second on xyz coordinates
		// Refer Geometry::setVertexSign() function for implementation.
	};
};


void Cell::setLargestFace(Face *face){
	this->largestFace = face;
};

Face* Cell::getLargestFace(){
	return this->largestFace;
};

int Cell::getNoOfBoundaryFaces(){
	return this->noOfBoundaryFaces;
};

void Cell::setNoOfBoundaryFaces(int NBF){
	this->noOfBoundaryFaces = NBF;
};

void Cell::addNoOfBoundaryFaces(int NBF){
	this->noOfBoundaryFaces += NBF;
};

void Cell::calculateCenter(){
	double sum_x=0.0,sum_y=0.0,sum_z=0.0;
	for(int i=0;i<this->getNoOfPoints();i++){
		sum_x += definingPoints[i]->getX();
		sum_y += definingPoints[i]->getY();
		sum_z += definingPoints[i]->getZ();
	}
	center.setX( sum_x/noOfPoints);
	center.setY( sum_y/noOfPoints);
	center.setZ( sum_z/noOfPoints);

	center.setId(this->getId());
}


Point* Cell::getCenter(){
	return &this->center;
};


void Cell::setReferenceCell(RefCell* refCell){
	this->refCell = refCell;
	this->refCellFlag = true;
};

RefCell* Cell::getReferenceCell(){
	assert (this->refCellFlag == true && "Cell::getReferenceCell() returns error. Make sure that nDOF number is assigned.\n");
	return this->refCell;
};

void Cell::setNDOF(unsigned int nDOF){
	this->nDOF = nDOF;
	this->nDOFFlag = true;
};

unsigned int Cell::getNDOF()const{
	assert (this->nDOFFlag == true && "Cell::getNDOF() returns error. Make sure that nDOF number is assigned.\n");
	return this->nDOF;
};

void Cell::setNQuad(unsigned int nQuad){
	this->nQuad = nQuad;
	this->nQuadFlag = true;
};

unsigned int Cell::getNQuad()const{
	assert (this->nQuadFlag == true && "Cell::getNQuad() returns error. Make sure that nQuad number is assigned.\n");
	return this->nQuad;
};


TensorO1<double>* Cell::getJacobian(){
	assert(this->JacobianFlag == true && "Cell::getJacobian() returns error. Make sure that the Jacobian is computed.\n");
	return &this->J;
};

TensorO3<double>* Cell::getInverseJacobian(){
	return &this->dRST_by_dXYZ;
};

void Cell::setVertexSign(int vertexNo, int coord, double sign){
	this->vertexSign.setValue(vertexNo, coord, sign);
};

double Cell::getVertexSign(int vertexNo, int coord){
	return this->vertexSign.getValue(vertexNo, coord);
};

TensorO2<double>* Cell::getVertexSignArray(){
	return &this->vertexSign;
};

TensorO2<double>* Cell::getQuadPointsGlobalLocation(){
	return &this->quadPointsGlobalLocation;
};

TensorO2<double>* Cell::getDOFPointsGlobalLocation(){
	return &this->DOFPointsGlobalLocation;
};

void Cell::assertVertexSigns(){
	enum coord{x,y,z};
	if (this->getCellType() == CellType::Hex){
		assert(this->getVertexSign(0,x)== -1.0 && this->getVertexSign(0,y)== -1.0 && this->getVertexSign(0,z)== -1.0);
		assert(this->getVertexSign(1,x)== -1.0 && this->getVertexSign(1,y)==  1.0 && this->getVertexSign(1,z)== -1.0);
		assert(this->getVertexSign(2,x)==  1.0 && this->getVertexSign(2,y)==  1.0 && this->getVertexSign(2,z)== -1.0);
		assert(this->getVertexSign(3,x)==  1.0 && this->getVertexSign(3,y)== -1.0 && this->getVertexSign(3,z)== -1.0);
		assert(this->getVertexSign(4,x)== -1.0 && this->getVertexSign(4,y)== -1.0 && this->getVertexSign(4,z)==  1.0);
		assert(this->getVertexSign(5,x)== -1.0 && this->getVertexSign(5,y)==  1.0 && this->getVertexSign(5,z)==  1.0);
		assert(this->getVertexSign(6,x)==  1.0 && this->getVertexSign(6,y)==  1.0 && this->getVertexSign(6,z)==  1.0);
		assert(this->getVertexSign(7,x)==  1.0 && this->getVertexSign(7,y)== -1.0 && this->getVertexSign(7,z)==  1.0);
	};

};

void Cell::mapRSTtoXYZ3DTensor(Point* localPt, Point* globalPt){
	// The shape functions are given in: http://textofvideo.nptel.ac.in/105106051/lec40.pdf, page number 8.

	// Check that the vertices have the correct sign 
	this->assertVertexSigns();

	// r,s,t are the coordinates in the local (cell) frame of reference, i.e. defined in the reference cell.
	double r = localPt->getX();
	double s = localPt->getY();
	double t = localPt->getZ();

	// Enum coord is used for assesing vertexSigns in x or y or z direction. Refer Cell::vertexSign and GeometryIO::assignVertexSigns
	// for details.

	enum coord {x,y,z};

	// XYZ is the local container for global coordinates
	double XYZ[3];
	for (Index i=0; i<3; i++){
		XYZ[i] = 0.0;
	};

	for(Index vertex=0; vertex<this->getNoOfPoints(); vertex++){

		double shape_function_i = 1.0/8.0 * (1.0 + this->getVertexSign(vertex, x) *r) * 
		                                    (1.0 + this->getVertexSign(vertex, y) *s) * 
						    (1.0 + this->getVertexSign(vertex, z) *t);

		double x_i = this->getDefiningPoint(vertex)->getX();
		double y_i = this->getDefiningPoint(vertex)->getY();
		double z_i = this->getDefiningPoint(vertex)->getZ();

		XYZ[0] += shape_function_i * x_i;
		XYZ[1] += shape_function_i * y_i;
		XYZ[2] += shape_function_i * z_i;
	};

	// Now the values are dumped in the global point 
	globalPt->setX(XYZ[0]);
	globalPt->setY(XYZ[1]);
	globalPt->setZ(XYZ[2]);
};



// Calculates Jacobian for hex cells 
// Takes r,s,t coordinates in reference cell coordinate system
// returns double value of the determinant of the Jacobian at r,s,t point
double Cell::calculateJacobian3DTensor(double r,double s,double t) {
	
	this->assertVertexSigns();

	double J[3][3];
	// set values to zero (initialization)
	for (Index i=0; i<3; i++){
		for (Index j=0; j<3; j++){
			J[i][j] = 0.0;
		};
	};

	double detJ;
	enum coord{x,y,z};
	
	for(Index vertex=0; vertex<this->getNoOfPoints(); vertex++){
		double signX = this->getVertexSign(vertex, x);
		double signY = this->getVertexSign(vertex, y);
		double signZ = this->getVertexSign(vertex, z);

		J[0][0] += signX * (1.0 + signY * s) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getX(); //dx_by_dr
		J[0][1] += signY * (1.0 + signX * r) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getX(); //dx_by_ds
		J[0][2] += signZ * (1.0 + signX * r) * (1.0 + signY * s) / 8.0 * this->getDefiningPoint(vertex)->getX(); //dx_by_dt

		J[1][0] += signX * (1.0 + signY * s) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getY(); //dy_by_dr
		J[1][1] += signY * (1.0 + signX * r) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getY(); //dy_by_ds
		J[1][2] += signZ * (1.0 + signX * r) * (1.0 + signY * s) / 8.0 * this->getDefiningPoint(vertex)->getY(); //dy_by_dt

		J[2][0] += signX * (1.0 + signY * s) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getZ(); //dz_by_dr
		J[2][1] += signY * (1.0 + signX * r) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getZ(); //dz_by_ds
		J[2][2] += signZ * (1.0 + signX * r) * (1.0 + signY * s) / 8.0 * this->getDefiningPoint(vertex)->getZ(); //dz_by_dt

	}
	
	// Direct formula for determinant		
	detJ = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
	return detJ;
};
	
		


// Calculates inverse Jacobian for the cell (not determinant)
// each term of the inverse Jacobian is a metric term for transformation
// The terms are stored in an array dRST_by_dXYZ. Dim: nQuad x 3 x 3 
// i.e. [ [dr/dx, dr/dy, dr/dz], 
//	  [ds/dx, ds/dy, ds/dz], 
//	  [dt/dx, dt/dy, dt/dz]]    matrix for nQuad points

// Takes r,s,t coordinates in reference cell coordinate system
// Fills the matrix for those coordinates (i.e. quadrature point)
double Cell::calculateInverseJacobianMatrix3DTensor(double r,double s,double t, Matrix<double> *dRST_by_dXYZ) {
	
	this->assertVertexSigns();
	assert(dRST_by_dXYZ->getSize0() == dRST_by_dXYZ->getSize1() and dRST_by_dXYZ->getSize1()==3);

	double J[3][3];
	// set values to zero (initialization)
	for (Index i=0; i<3; i++){
		for (Index j=0; j<3; j++){
			J[i][j] = 0.0;
		};
	};

	double detJ;
	enum coord{x,y,z};
	
	for(Index vertex=0; vertex<this->getNoOfPoints(); vertex++){
		double signX = this->getVertexSign(vertex, x);
		double signY = this->getVertexSign(vertex, y);
		double signZ = this->getVertexSign(vertex, z);

		J[0][0] += signX * (1.0 + signY * s) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getX(); //dx_by_dr
		J[0][1] += signY * (1.0 + signX * r) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getX(); //dx_by_ds
		J[0][2] += signZ * (1.0 + signX * r) * (1.0 + signY * s) / 8.0 * this->getDefiningPoint(vertex)->getX(); //dx_by_dt
                                                                                                                                   
		J[1][0] += signX * (1.0 + signY * s) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getY(); //dy_by_dr
		J[1][1] += signY * (1.0 + signX * r) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getY(); //dy_by_ds
		J[1][2] += signZ * (1.0 + signX * r) * (1.0 + signY * s) / 8.0 * this->getDefiningPoint(vertex)->getY(); //dy_by_dt
                                                                                                                                   
		J[2][0] += signX * (1.0 + signY * s) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getZ(); //dz_by_dr
		J[2][1] += signY * (1.0 + signX * r) * (1.0 + signZ * t) / 8.0 * this->getDefiningPoint(vertex)->getZ(); //dz_by_ds
		J[2][2] += signZ * (1.0 + signX * r) * (1.0 + signY * s) / 8.0 * this->getDefiningPoint(vertex)->getZ(); //dz_by_dt

	}
	
	// Direct formula for determinant		
	detJ = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);

	// Direct formula for the inverse of J 
	// Refer: https://www.dr-lex.be/random/matrix-inv.html
	
	// Computing each of the Metric terms
	double drBydx = 1.0/detJ*( J[2][2]*J[1][1] - J[2][1]*J[1][2]) ;
	double drBydy = -1.0/detJ*(J[2][2]*J[0][1] - J[2][1]*J[0][2]) ;
	double drBydz = 1.0/detJ*( J[1][2]*J[0][1] - J[1][1]*J[0][2])  ;

	double dsBydx = -1.0/detJ*(J[2][2]*J[1][0] - J[2][0]*J[1][2])  ;
	double dsBydy = 1.0/detJ*( J[2][2]*J[0][0] - J[2][0]*J[0][2])  ;
	double dsBydz = -1.0/detJ*(J[1][2]*J[0][0] - J[1][0]*J[0][2])  ;

	double dtBydx = 1.0/detJ*( J[2][1]*J[1][0] - J[2][0]*J[1][1])  ;
	double dtBydy = -1.0/detJ*(J[2][1]*J[0][0] - J[2][0]*J[0][1])  ;
	double dtBydz = 1.0/detJ*( J[1][1]*J[0][0] - J[1][0]*J[0][1])  ;


	// The inverse Jacobian (i.e. metric terms for derivatives given as follows)

	enum IJlocal  {R, S, T}; 
	enum IJglobal {X, Y, Z}; 

	dRST_by_dXYZ->setValue(R, X, drBydx);
	dRST_by_dXYZ->setValue(R, Y, drBydy);
	dRST_by_dXYZ->setValue(R, Z, drBydz);

	dRST_by_dXYZ->setValue(S, X, dsBydx);
	dRST_by_dXYZ->setValue(S, Y, dsBydy);
	dRST_by_dXYZ->setValue(S, Z, dsBydz);

	dRST_by_dXYZ->setValue(T, X, dtBydx);
	dRST_by_dXYZ->setValue(T, Y, dtBydy);
	dRST_by_dXYZ->setValue(T, Z, dtBydz);
};



Matrix<double>* Cell::getMassMatrix(){
	return &this->M;
};

Matrix<double>* Cell::getInverseMassMatrix(){
	return &this->MI;
};

Matrix<double>* Cell::getVandermondeMatrix(){
	return &this->V;
};

Matrix<double>* Cell::getDxMatrix(){
	return &this->Dx;
};

Matrix<double>* Cell::getDyMatrix(){
	return &this->Dy;
};

Matrix<double>* Cell::getDzMatrix(){
	return &this->Dz;
};

Matrix<double>* Cell::getDtildexMatrix(){
	return &this->Dtildex;
};

Matrix<double>* Cell::getDtildeyMatrix(){
	return &this->Dtildey;
};

Matrix<double>* Cell::getDtildezMatrix(){
	return &this->Dtildez;
};

TensorO3<double>* Cell::getFMatrix(){
	return &this->F;
};

TensorO3<double>* Cell::getFtildeMatrix(){
	return &this->Ftilde;
};

Matrix<double>* Cell::getFilterMatrix(){
	return &this->Filter;
};

void Cell::resizeMatrices(unsigned int Size){
	assert(nDOFFlag == true and nQuadFlag == true and this->getNoOfFaces()>0); 

	// Re-size all the matrices and set all entries to zero 

	// Mass matrix
	// Resize the matrix
	this->getMassMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getMassMatrix()->getSize0() == this->getMassMatrix()->getSize1() and this->getMassMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getMassMatrix()->setAll(0.0);
	

	// Inverse Mass matrix
	this->getInverseMassMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getInverseMassMatrix()->getSize0() == this->getInverseMassMatrix()->getSize1() and this->getInverseMassMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getInverseMassMatrix()->setAll(0.0);
	
	// Vandermonde matrix
	this->getVandermondeMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getVandermondeMatrix()->getSize0() == this->getVandermondeMatrix()->getSize1() and this->getVandermondeMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getVandermondeMatrix()->setAll(0.0);

	// Dx matrix
	this->getDxMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getDxMatrix()->getSize0() == this->getDxMatrix()->getSize1() and this->getDxMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getDxMatrix()->setAll(0.0);


	// Dy matrix
	this->getDyMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getDyMatrix()->getSize0() == this->getDyMatrix()->getSize1() and this->getDyMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getDyMatrix()->setAll(0.0);


	// Dz matrix
	this->getDzMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getDzMatrix()->getSize0() == this->getDzMatrix()->getSize1() and this->getDzMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getDzMatrix()->setAll(0.0);


	// DTildeX matrix
	this->getDtildexMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getDtildexMatrix()->getSize0() == this->getDtildexMatrix()->getSize1() and this->getDtildexMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getDtildexMatrix()->setAll(0.0);

	// DTildeY matrix
	this->getDtildeyMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getDtildeyMatrix()->getSize0() == this->getDtildeyMatrix()->getSize1() and this->getDtildeyMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getDtildeyMatrix()->setAll(0.0);

	// DTildeZ matrix
	this->getDtildezMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getDtildezMatrix()->getSize0() == this->getDtildezMatrix()->getSize1() and this->getDtildezMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getDtildezMatrix()->setAll(0.0);

	// Flux matrix
	this->getFMatrix()->setSize(this->getNoOfFaces(), Size, Size);
	// Confirm the size
	assert(this->getFMatrix()->getSize0() == this->getNoOfFaces() and  this->getFMatrix()->getSize1() == this->getFMatrix()->getSize2() and this->getFMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getFMatrix()->setAll(0.0);

	// FTilde matrix
	this->getFtildeMatrix()->setSize(this->getNoOfFaces(), Size, Size);
	// Confirm the size
	assert(this->getFtildeMatrix()->getSize0() == this->getNoOfFaces() and  this->getFtildeMatrix()->getSize1() == this->getFtildeMatrix()->getSize2() and this->getFtildeMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getFtildeMatrix()->setAll(0.0);

	// Filter matrix
	this->getFilterMatrix()->setSize(Size);
	// Confirm the size
	assert(this->getFilterMatrix()->getSize0() == this->getFilterMatrix()->getSize1() and this->getFilterMatrix()->getSize1() == this->getNDOF());
	// Initialize with zero
	this->getFilterMatrix()->setAll(0.0);

};



double Cell::getVariable(int rkStep, int varNo, int DOF){
	return this->variable.getValue(rkStep, varNo, DOF);
};

TensorO3<double>* Cell::getVariableArray(){
	return &this->variable;
};

TensorO2<double>* Cell::getVariableRefArray(){
	return &this->variableRef;
};
TensorO3<double>* Cell::getVariableResidualArray(){
	return &this->variableResidual;
};

TensorO2<double>* Cell::getFluxVectorx(){
	return &this->fluxVectorx;
};

TensorO2<double>* Cell::getFluxVectory(){
	return &this->fluxVectory;
};

TensorO2<double>* Cell::getFluxVectorz(){
	return &this->fluxVectorz;
};

TensorO3<double>* Cell::getFluxStarVector(){
	return &this->fluxStarVector;
};


/* Printing Functions */


//Mass matrix

void Cell::printMassMatrix(){
	cout <<"\nM: \n";
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << M.getValue(i,j) << " ";
		};
		cout << endl;
	};
};

void Cell::printMassMatrix(ofstream &file){
	file <<"\nM: \n";
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << M.getValue(i,j) << " ";
		};
		file << endl;
	};
};


void Cell::printInverseMassMatrix(){
	cout <<"\nM_inv: \n";
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << MI.getValue(i, j) << " ";
		};
		cout << endl;
	};
};


void Cell::printInverseMassMatrix(ofstream &file){
	file <<"\nM_inv: \n";
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << MI.getValue(i, j) << " ";
		};
		file << endl;
	};
};


//Vandermonde matrix
void Cell::printVandermondeMatrix(){
	cout <<"\nVandermonde: \n";
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << V.getValue(i,j) << " ";
		};
		cout << endl;
	};
};

void Cell::printVandermondeMatrix(ofstream &file){
	file <<"\nVandermonde: \n";
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << V.getValue(i,j) << " ";
		};
		file << endl;
	};
};


//Differential matrix
void Cell::printDifferentialMatrix(){
	cout << "\nDx: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << Dx.getValue(i,j) << " ";
		};
		cout << endl;
	};
	cout << "\nDy: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << Dy.getValue(i,j) << " ";
		};
		cout << endl;
	};
	cout << "\nDz: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << Dz.getValue(i,j) << " ";
		};
		cout << endl;
	};
};


//Differential matrix
void Cell::printDifferentialMatrix(ofstream &file){
	file << "\nDx: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << Dx.getValue(i,j) << " ";
		};
		file << endl;
	};
	file << "\nDy: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << Dy.getValue(i,j) << " ";
		};
		file << endl;
	};
	file << "\nDz: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << Dz.getValue(i,j) << " ";
		};
		file << endl;
	};
};



// Differential matrix multiplied with M_inv
void Cell::printDTildeMatrix(){
	cout << "\nDTildeX: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << Dtildex.getValue(i,j) << " ";
		};
		cout << endl;
	};
	cout << "\nDTildeY: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << Dtildey.getValue(i,j) << " ";
		};
		cout << endl;
	};
	cout << "\nDTildeZ: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << Dtildez.getValue(i,j) << " ";
		};
		cout << endl;
	};
};

// Differential matrix multiplied with M_inv
void Cell::printDTildeMatrix(ofstream &file){
	file << "\nDTildeX: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << Dtildex.getValue(i,j) << " ";
		};
		file << endl;
	};
	file << "\nDTildeY: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << Dtildey.getValue(i,j) << " ";
		};
		file << endl;
	};
	file << "\nDTildeZ: " << endl;
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << Dtildez.getValue(i,j) << " ";
		};
		file << endl;
	};
};



//Flux matrix
void Cell::printFluxMatrix(){
	cout << "\nFlux matrix:\n";

	for (int f=0; f<noOfFaces; f++){
		cout << "Face: " << f << endl;
		for (int i=0; i<nDOF; i++){
			for (int j=0; j<nDOF; j++){
				cout << F.getValue(f,i,j) << " ";
			};
			cout << endl;
		};
		cout << endl;
	};
};

void Cell::printFluxMatrix(ofstream &file){
	file << "\nFlux matrix:\n";

	for (int f=0; f<noOfFaces; f++){
		file << "Face: " << f << endl;
		for (int i=0; i<nDOF; i++){
			for (int j=0; j<nDOF; j++){
				file << F.getValue(f,i,j) << " ";
			};
			file << endl;
		};
		file << endl;
	};
};

//Flux matrix
void Cell::printFTildeMatrix(){
	cout << "\nF-tilde (M_inv * F) matrix:\n";

	for (int f=0; f<noOfFaces; f++){
		cout << "Face: " << f << endl;
		for (int i=0; i<nDOF; i++){
			for (int j=0; j<nDOF; j++){
				cout << Ftilde.getValue(f,i,j) << " ";
			};
			cout << endl;
		};
		cout << endl;
	};
};

//Flux matrix
void Cell::printFTildeMatrix(ofstream &file){
	file << "\nF-tilde (M_inv * F) matrix:\n";

	for (int f=0; f<noOfFaces; f++){
		file << "Face: " << f << endl;
		for (int i=0; i<nDOF; i++){
			for (int j=0; j<nDOF; j++){
				file << Ftilde.getValue(f,i,j) << " ";
			};
			file << endl;
		};
		file << endl;
	};
};

// Filter matrix
void Cell::printFilterMatrix(){
	cout <<"\nFilter Matrix: \n";
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			cout << Filter.getValue(i,j) << " ";
		};
		cout << endl;
	};
};

void Cell::printFilterMatrix(ofstream &file){
	file <<"\nFilter Matrix: \n";
	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			file << Filter.getValue(i,j) << " ";
		};
		file << endl;
	};
};









/*******************************************TESTS**************************************************/


namespace Test{

	void TestCell(){
		cout << "\nRunning test on class Cell... " << endl;

		TestCellInitialization();
		TestDefiningPointsForCell();
		TestDefiningFacesForCell();

		cout << "Test::TestCell() passed.\n";
	};

	void TestCellInitialization(){
		// Run the test on the dataset 1000 times using random values
		for (Index i=0; i<1000; i++){
			int random = getRandom(-10000, 10000);
			Cell cell(random);
			int grandom = cell.getId();
			assert(random == grandom && "Test::TestCellInitialization() returns error.\n");
		};
	};

	void TestDefiningPointsForCell(){
		// Run the test on the dataset 1000 times using random values
		for (Index i=0; i<1000; i++){
			Point points[8];
			int IDP[8], IDF[8];
			double XP[8], XF[8];
			double YP[8], YF[8];
			double ZP[8], ZF[8];

			int random = getRandom(-10000, 10000);
			Cell cell(random);

			int noOfPts = getRandom(4,8);
			if (noOfPts == 7){
				noOfPts = 8;
			};

			for (Index i=0; i<noOfPts; i++){
				double x1, y1, z1; int id1;
				x1 = getRandom(-1000, 1000);
				y1 = getRandom(-1000, 1000);
				z1 = getRandom(-1000, 1000);
				id1 = getRandom(-1000, 1000);
				points[i].setId(id1);
				points[i].setX(x1);
				points[i].setY(y1);
				points[i].setZ(z1);

				IDP[i] = id1;
				XP[i] = x1;
				YP[i] = y1;
				ZP[i] = z1;

				cell.addDefiningPoint(&points[i]);
			};

			int Np = cell.getNoOfPoints();
			assert(Np == noOfPts && "Test::TestDefiningPoints() in Cell.cpp returns error. Cell::getNoOfPoints() returns wrong value.\n");
			for (Index i=0; i<cell.getNoOfPoints(); i++){
				IDF[i] = cell.getDefiningPoint(i)->getId();

				XF[i] = cell.getDefiningPoint(i)->getX();
				YF[i] = cell.getDefiningPoint(i)->getY();
				ZF[i] = cell.getDefiningPoint(i)->getZ();

				assert(IDP[i] == IDF[i] and "Test::TestDefiningPoints() in Cell.cpp returns error. Cell::getDefiningPoint() does not return the correct id.\n");
				assert(XP[i] == XF[i] and YP[i] == YF[i] and ZP[i] == ZF[i] && "Test::TestDefiningPoints() in Cell.cpp returns error. Cell::getDefiningPoint() does not return the correct xyz coordinate.\n");
			}



			// test the findPoint() method
			bool flag = false;
			flag = cell.findPoint(&points[getRandom(1,Np)]);

			assert(flag == true && "Test::TestDefiningPoints() in cell.cpp returns error. Function Cell::findPoint() does not work properly.\n");

			Point xyz(100);
			flag = cell.findPoint(&xyz);

			assert(flag == false && "Test::TestDefiningPoints() in cell.cpp returns error. Function Cell::findPoint() does not work properly.\n");
		};

	};

	void TestDefiningFacesForCell(){
		Cell cell(getRandom(100,200));
		int nfa = getRandom(4,6);
		cell.setNoOfFaces(nfa);
		int nf = cell.getNoOfFaces();

		assert(nf == nfa);


		for (Index i=0; i<1000; i++){
			Cell cell(i);

			int nfa = getRandom(4,6);

			Face *faces[nfa];

			for (Index j=0; j<nfa; j++){
				Face face(j);
				faces[j] = &face;
				cell.addDefiningFace(faces[j]);
			};


			for (Index j=0; j<nfa; j++){
				assert(faces[j]->getId() == cell.getDefiningFace(j)->getId());
			};

			// At this time, the orderedFaces array just contains numbers 0 to 8 arranged sequencially
			// The next routine rearranges the faces based on the correct order

		};

		// Test for finding defining points and defining Faces

		for (Index i=0; i<1000; i++){
			Cell cell(i);

			// Set 8 points with ids 0 to 7
			Point *points;
			points = new Point[8];

			for (Index j=0; j<8; j++){
				points[j].setId(j);
			};

			// Set 6 faces with ids 0 to 5
			Face *faces;

			faces = new Face[6];

			for (Index j=0; j<6; j++){
				faces[j].setId(j);
			};


			// Zeroth face gets the first four points in that sequence
			for (Index j=0; j<4; j++){
				faces[0].addDefiningPoint(&points[j]);
			};
			faces[0].setFaceType();

			// Now create an array of numbers 0 to 5 and shuffle it (100 times)
			TensorO1<int> order(5);
			for (Index j=0; j<5; j++){
				order.setValue(j,j+1);
			};
			order.randomize();	

			// Now randomly arrange face ids 

			// Face sharing points 0 and 1 with the base face 
			faces[order.getValue(0)].addDefiningPoint(&points[0]);
			faces[order.getValue(0)].addDefiningPoint(&points[4]);
			faces[order.getValue(0)].addDefiningPoint(&points[5]);
			faces[order.getValue(0)].addDefiningPoint(&points[1]);
			faces[order.getValue(0)].setFaceType();

			// Face sharing points 1 and 2 with the base face 
			faces[order.getValue(1)].addDefiningPoint(&points[1]);
			faces[order.getValue(1)].addDefiningPoint(&points[5]);
			faces[order.getValue(1)].addDefiningPoint(&points[6]);
			faces[order.getValue(1)].addDefiningPoint(&points[2]);
			faces[order.getValue(1)].setFaceType();

			// Face sharing points 2 and 3 with the base face 
			faces[order.getValue(2)].addDefiningPoint(&points[2]);
			faces[order.getValue(2)].addDefiningPoint(&points[6]);
			faces[order.getValue(2)].addDefiningPoint(&points[7]);
			faces[order.getValue(2)].addDefiningPoint(&points[3]);
			faces[order.getValue(2)].setFaceType();

			// Face sharing points 3 and 0 with the base face 
			faces[order.getValue(3)].addDefiningPoint(&points[3]);
			faces[order.getValue(3)].addDefiningPoint(&points[7]);
			faces[order.getValue(3)].addDefiningPoint(&points[4]);
			faces[order.getValue(3)].addDefiningPoint(&points[0]);
			faces[order.getValue(3)].setFaceType();

			// Face not sharing any points with the base face (i.e. top face)
			faces[order.getValue(4)].addDefiningPoint(&points[4]);
			faces[order.getValue(4)].addDefiningPoint(&points[7]);
			faces[order.getValue(4)].addDefiningPoint(&points[6]);
			faces[order.getValue(4)].addDefiningPoint(&points[5]);
			faces[order.getValue(4)].setFaceType();


			// Now add these faces in to the cell a random order

			TensorO1<int> order2(5);
			for (Index j=0; j<5; j++){
				order2.setValue(j,order.getValue(j));
			};
			order2.randomize();

			cell.addDefiningFace(&faces[0]);
			cell.addDefiningFace(&faces[order2.getValue(0)]);
			cell.addDefiningFace(&faces[order2.getValue(1)]);
			cell.addDefiningFace(&faces[order2.getValue(2)]);
			cell.addDefiningFace(&faces[order2.getValue(3)]);
			cell.addDefiningFace(&faces[order2.getValue(4)]);


			// After all the faces are arranged, the cell is finally constructed.
			cell.findDefiningPoints();
			assert(cell.getNoOfPoints() == 8);

			cell.setCellType(); 		//implicitly finds the order of points and faces
			assert(cell.getCellType() == CellType::Hex);

			// The final sequence of the faces should be:
			// Base: 0,  (containing points 0 1 2 3)
			// Top : 1,  (containing points 4 5 6 7), with 4 opposite of 0, 5 opposite of 1 etc.
			// Front:2,  (sharing pts 0 and 1 with base) 
			// LHS:  3,  (sharing pts 1 and 2 with base)
			// Back: 4,  (sharing pts 2 and 3 with base)
			// RHS:  5,  (sharing pts 3 and 4 with base)

			assert(cell.getDefiningFace(0)->getId() == 0);
			assert(cell.getDefiningFace(1)->getId() == order.getValue(4));
			assert(cell.getDefiningFace(2)->getId() == order.getValue(0));
			assert(cell.getDefiningFace(3)->getId() == order.getValue(1));
			assert(cell.getDefiningFace(4)->getId() == order.getValue(2));
			assert(cell.getDefiningFace(5)->getId() == order.getValue(3));


			// Next, make sure that the points are returned in the correct sequence:
			for (Index j=0; j<8; j++){
				assert(cell.getDefiningPoint(j)->getId() == j);
			};

			delete[] faces;
			delete[] points;
		};
	};
};

