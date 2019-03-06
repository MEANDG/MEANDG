#include "Face.h"


// Constructor definition 
Face::Face(int Id){
	this->id = Id;
	this->noOfPoints = 0;
	this->area = EPS;
	
	this->owner = NULL;
	this->neighbour = NULL;

	// Flags
	this->faceTypeFlag = false;		// Face type is not assigned during initialization.
	this->boundaryFlag = true;		// Each face is assumed to be a boundary face until neighbour is found out.
	this->normalFlag = false;		// NormalFlag is set to true after normal and tangents are computed.
	this->nDOFFlag  = false;		// nDOFFlag is set to true after assigning number of DOF locations
	this->nQuadFlag  = false;		// nQuadFlag is set to true after assigning the number of quadrature points


	// Default order of the defining points
	for (Index i=0; i<4; i++){
		this->orderedPoints[i] = i;
	};

	this->refFaceFlag = false;
	this->JacobianFlag = false;

	this->getRotationMatrixParallelToXY()->setSize(3);
	this->getInverseRotationMatrixParallelToXY()->setSize(3);
	this->getRotationMatrixNormal()->setSize(3);
	this->getInverseRotationMatrixNormal()->setSize(3);
}


Face::~Face(){
};


void Face::setId(unsigned int Id){
	this->id = Id;
};

int Face::getId(){
	return this->id;
};

void Face::setName(string name){
	this->name = name;
};

string Face::getName(){
	return this->name;
};

unsigned int Face::getNoOfPoints() const{
	assert(this->noOfPoints >= 3 && "Face::getNoOfPoints returns error. Make sure that all the points are added\n.");
	return this->noOfPoints;
};

Point* Face::getDefiningPoint(int i) const{
	assert(i < this->noOfPoints && "Face::getDefiningPoint() returns error. Argument index greater than total number of defining points.\n");
	return this->definingPoints[this->orderedPoints[i]];
};


void Face::addDefiningPoint(Point *definingPoint){
	assert(this->noOfPoints <= 3 && "Face::addDefiningPoint() returns error. Face already has atleast 4 points.\n");
	assert(this->faceTypeFlag == false && "Face::addDefiningPoint() returns error. The face type is already finalized.\n");
	this->definingPoints[this->noOfPoints] = definingPoint;
	this->noOfPoints ++;
};


void Face::setFaceType(){
	assert(this->faceTypeFlag == false && "Face::setFaceType() returns error. The face type is already finalized.\n");
	assert(this->noOfPoints == 3 or this->noOfPoints == 4 && "Face::setFaceType() returns error. Build the face completely first.\n");

	if (this->noOfPoints == 3){
		facetype = FaceType::Tri;
	}
	else{
		facetype = FaceType::Quad;
	};
	this->faceTypeFlag = true; 	//once this flag is set to true, face type can not be reset or new point can not be added.
};


FaceType::faceType Face::getFaceType() const{
	assert(this->faceTypeFlag == true && "Face::setFaceType() returns error. Set the facetype first.\n");
	return this->facetype;
};


void Face::setOwnerCell(Cell *cell){
	assert(this->owner == NULL && "Face::setOwnerCell() returns error. Owner cell already assigned.\n");	
	this->owner = cell;
};

Cell* Face::getOwnerCell(){
	return this->owner;
};

void Face::setNeighbourCell(Cell *cell){
	assert(this->neighbour == NULL && "Face::setNeighbourCell() returns error. Neighbour cell already assigned.\n");	
	this->neighbour = cell;
	this->boundaryFlag = false;
};

Cell* Face::getNeighbourCell(){
	return this->neighbour;
};

bool Face::isBoundary(){
	assert(this->owner != NULL && "Face::isBoundary() returns error. Make sure that both the owner and the neighbour files are processed.\n");
	return this->boundaryFlag;
};


void Face::calculateNormal(){
	assert(this->faceTypeFlag == true && "Face::calculateNormal() returns error. Make sure that face is complete and faceTypeFlag is set to true.\n");

	/// Assign space for normal and tangent tensors
	this->normal.setSize(3);
	this->tangent1.setSize(3);
	this->tangent2.setSize(3);

	double x[3],y[3],z[3],a[3],b[3];
	double Normal[3], Tangent1[3], Tangent2[3];
	
	for(int i=0;i<3;i++){
		x[i] = definingPoints[i]->getX();
		y[i] = definingPoints[i]->getY();
		z[i] = definingPoints[i]->getZ();
	}
	a[0] = x[1]-x[0];	
	a[1] = y[1]-y[0];	
	a[2] = z[1]-z[0];	

	b[0] = x[2]-x[0];	
	b[1] = y[2]-y[0];	
	b[2] = z[2]-z[0];	
	
	Normal[0] = a[1]*b[2]-a[2]*b[1];
	Normal[1] = a[2]*b[0]-a[0]*b[2];
	Normal[2] = a[0]*b[1]-a[1]*b[0];
	
	double Normal_mag = pow((pow(Normal[0],2)+pow(Normal[1],2)+pow(Normal[2],2)),0.5);
	Normal[0] /= Normal_mag;
	Normal[1] /= Normal_mag;
	Normal[2] /= Normal_mag;
	
	// set normal value
	for (Index i =0; i<3; i++){
		this->normal.setValue(i, Normal[i]);
	};

	// computing Tangent 1

	Tangent1[0] = Normal[1];
	Tangent1[1] = -Normal[0];
	Tangent1[2] = 0.0;

	if (abs(Normal[0]) < NANO and abs(Normal[1]) < NANO){
		Tangent1[0] = 0.0;
		Tangent1[1] = Normal[2];
		Tangent1[2] = -Normal[1];
	};

	// Normalize Tangent 1 to get a unit vector
	double tanMag;
	tanMag = pow((pow(Tangent1[0],2)+pow(Tangent1[1],2)+pow(Tangent1[2],2)),0.5);
	Tangent1[0] /= tanMag;
	Tangent1[1] /= tanMag;
	Tangent1[2] /= tanMag;

	// set tangent 1 value
	for (Index i = 0; i<3; i++){
		this->tangent1.setValue(i, Tangent1[i]);
	};


	// take cross product of Normal and Tangent1 to get the mutually perpendicular Tangent2
	Tangent2[0] = Normal[1]*Tangent1[2] - Normal[2]*Tangent1[1];
	Tangent2[1] = -(Normal[0]*Tangent1[2] - Normal[2]*Tangent1[0]);
	Tangent2[2] = Normal[0]*Tangent1[1] - Normal[1]*Tangent1[0];

	// Normalize Tangent 2 to get a unit vector
	tanMag = pow((pow(Tangent2[0],2)+pow(Tangent2[1],2)+pow(Tangent2[2],2)),0.5);
	Tangent2[0] /= tanMag;
	Tangent2[1] /= tanMag;
	Tangent2[2] /= tanMag;

	// set tangent 2 value
	for (Index i = 0; i<3; i++){
		this->tangent2.setValue(i, Tangent2[i]);
	};


	// make sure that the Tangents are computed correctly
	assert(abs(Math::dot(normal,tangent1)) < NANO && "Face::calculateNormal() returns error. Check the Tangent vector computation in Face.cpp");
	assert(abs(Math::dot(normal,tangent2)) < NANO && "Face::calculateNormal() returns error. Check the Tangent vector computation in Face.cpp");
	assert(abs(Math::dot(tangent1,tangent2)) < NANO &&"Face::calculateNormal() returns error. Check the Tangent vector computation in Face.cpp");

	// Normal flag is set to true
	this->normalFlag = true;
	
	// Order the defining points based on the normal orientation
	this->orderDefiningPoints();

	this->vertexSign.setSize(this->getNoOfPoints(), 2);	// first index on defining points, second on xy coordinates

};


TensorO1<double>* Face::getNormal(){
	assert(this->normalFlag == true && "Face::getNormal() returns error. Make sure that the normal is computed.\n");
	return &this->normal;
};

TensorO1<double>* Face::getTangent1(){
	assert(this->normalFlag == true && "Face::getTangent1() returns error. Make sure that the tangent1 is computed.\n");
	return &this->tangent1;
};

TensorO1<double>* Face::getTangent2(){
	assert(this->normalFlag == true && "Face::getTangent2() returns error. Make sure that the tangent2 is computed.\n");
	return &this->tangent2;
};


/// Order the defining points in the correct order. The order is maintained such that the normal follows the right hand thumb rule.
/// i.e. if the fingers are curled in the sequence of the defining points, then the thumb points in the direction of the normal.

void Face::orderDefiningPoints(){
	assert(this->normalFlag == true && "Face::orderDefiningPoints() returns error. Make sure that the normal vector is computed.\n");

	if (this->facetype == FaceType::Tri){
		double p0[3], p1[3], p2[3];

		p0[0] = this->definingPoints[0]->getX();
		p0[1] = this->definingPoints[0]->getY();
		p0[2] = this->definingPoints[0]->getZ();

		p1[0] = this->definingPoints[1]->getX();
		p1[1] = this->definingPoints[1]->getY();
		p1[2] = this->definingPoints[1]->getZ();

		p2[0] = this->definingPoints[2]->getX();
		p2[1] = this->definingPoints[2]->getY();
		p2[2] = this->definingPoints[2]->getZ();

		// Now, find out the vectors originating at 0 and going to points 1, 2 and 3
		TensorO1<double> V01(3), V02(3);

		for (Index i=0; i<3; i++){
			V01.setValue(i, p1[i] - p0[i]);
			V02.setValue(i, p2[i] - p0[i]);
		};

		TensorO1<double> normalVector(3);
		Math::cross(V01, V02, &normalVector);

		if (Math::dot(normalVector, this->normal) > 0){
			this->orderedPoints[0] = 0;
			this->orderedPoints[1] = 1;
			this->orderedPoints[2] = 2;
		}
		else{
			this->orderedPoints[0] = 0;
			this->orderedPoints[1] = 2;
			this->orderedPoints[2] = 1;
		};
	}
	else if (this->facetype == FaceType::Quad){
		// Assuming vertices are named as 0, 1, 2 and 3 in any order
		// The required order of the face vertices is:
		/*
		       3*_______________________*2
                        |                       |
                        |                       |
                        |                       |
                        |                       |
                        |                       |
                        |                       |
                        |                       |
		   	*_______________________*
			0			1

			      i.e. right hand thumb rule. Normal coming out of the screen.
		*/

		double p0[3], p1[3], p2[3], p3[3];

		p0[0] = this->definingPoints[0]->getX();
		p0[1] = this->definingPoints[0]->getY();
		p0[2] = this->definingPoints[0]->getZ();

		p1[0] = this->definingPoints[1]->getX();
		p1[1] = this->definingPoints[1]->getY();
		p1[2] = this->definingPoints[1]->getZ();

		p2[0] = this->definingPoints[2]->getX();
		p2[1] = this->definingPoints[2]->getY();
		p2[2] = this->definingPoints[2]->getZ();

		p3[0] = this->definingPoints[3]->getX();
		p3[1] = this->definingPoints[3]->getY();
		p3[2] = this->definingPoints[3]->getZ();

		// Now, find out the vectors originating at 0 and going to points 1, 2 and 3
		TensorO1<double> V01(3), V02(3), V03(3);

		for (Index i=0; i<3; i++){
			V01.setValue(i, p1[i] - p0[i]);
			V02.setValue(i, p2[i] - p0[i]);
			V03.setValue(i, p3[i] - p0[i]);
		};

		// Next, find angles between three pairs formed by three vectors (since 0 base is common)
		double theta102 /*angle between V01 and V02*/, 
		       theta203 /*angle between V01 and V02*/, 
		       theta103 /*angle between V01 and V02*/;

		// Note: absolute value of the angle is considered sicne sign doesn't matter
		theta102 = abs(acos(Math::dot(V01,V02) / (V01.getL2Norm() * V02.getL2Norm())));
		theta103 = abs(acos(Math::dot(V01,V03) / (V01.getL2Norm() * V03.getL2Norm())));
		theta203 = abs(acos(Math::dot(V02,V03) / (V02.getL2Norm() * V03.getL2Norm())));

		
		double pa1[3], pa2[3], po[3]; 	//pa1: padjacent1,   pa2: padjacent2,    po: popposite

		double thetaMax = max(max(theta102, theta103), theta203);

		if ( thetaMax == theta102){
			// p1 and p2 are adjacent to p0; p3 is opposite

			TensorO1<double> normalVector(3);
			Math::cross(V01, V02, &normalVector);

			if (Math::dot(normalVector, this->normal) > 0){
				this->orderedPoints[0] = 0;
				this->orderedPoints[1] = 1;
				this->orderedPoints[2] = 3;
				this->orderedPoints[3] = 2;
			}
			else{
				this->orderedPoints[0] = 0;
				this->orderedPoints[1] = 2;
				this->orderedPoints[2] = 3;
				this->orderedPoints[3] = 1;
			};

		}
		else if (thetaMax == theta103){
			// p1 and p3 are adjacent to p0; p2 is opposite

			TensorO1<double> normalVector(3);
			Math::cross(V01, V03, &normalVector);

			if (Math::dot(normalVector, this->normal) > 0){
				this->orderedPoints[0] = 0;
				this->orderedPoints[1] = 1;
				this->orderedPoints[2] = 2;
				this->orderedPoints[3] = 3;
			}
			else{
				this->orderedPoints[0] = 0;
				this->orderedPoints[1] = 3;
				this->orderedPoints[2] = 2;
				this->orderedPoints[3] = 1;
			};
		}
		else if (thetaMax == theta203){
			// p2 and p3 are adjacent to p0; p1 is opposite

			TensorO1<double> normalVector(3);
			Math::cross(V02, V03, &normalVector);

			if (Math::dot(normalVector, this->normal) > 0){
				this->orderedPoints[0] = 0;
				this->orderedPoints[1] = 2;
				this->orderedPoints[2] = 1;
				this->orderedPoints[3] = 3;
			}
			else{
				this->orderedPoints[0] = 0;
				this->orderedPoints[1] = 3;
				this->orderedPoints[2] = 1;
				this->orderedPoints[3] = 2;
			};
		};

	};
};

int Face::getOrderedPoints(int index){
	return this->orderedPoints[index];	
};


void Face::calculateArea(){
	assert(this->faceTypeFlag == true && "Face::calculateArea() returns error. Make sure that all the points are added and face type is updated.\n");
	// First make sure that the points are ordered correctly. (orderDefiningPoints() is automatically called from calculateNormal())
	assert(this->normalFlag == true && "Face::calculateArea() returns error. Make sure that the normal vector is computed.\n");

	// Now the area computation.
	this->area = 0.0;

	if (this->facetype == FaceType::Tri){
		this->area = Math::calculateAreaOfTriangle(this->getDefiningPoint(0), this->getDefiningPoint(1), this->getDefiningPoint(2));
	}
	else if (this->facetype == FaceType::Quad){

		double Area1 = Math::calculateAreaOfTriangle(this->getDefiningPoint(0), this->getDefiningPoint(1), this->getDefiningPoint(2));
		double Area2 = Math::calculateAreaOfTriangle(this->getDefiningPoint(2), this->getDefiningPoint(3), this->getDefiningPoint(0));

		this->area = Area1 + Area2;
	};
	
};

double Face::getArea()const{
	assert(this->area > FEMTO && "Face::getArea() returns error. Please make sure that the area is computed.\n");
	return this->area;
};



void Face::calculateCenter(){
	assert(this->faceTypeFlag == true && "Face::setFaceType() returns error. Set the facetype first.\n");

	double sum_x=0.0,sum_y=0.0,sum_z=0.0;
	for(int i=0;i<this->noOfPoints;i++){
		sum_x += definingPoints[i]->getX();
		sum_y += definingPoints[i]->getY();
		sum_z += definingPoints[i]->getZ();
	}

	center.setId(this->id);
	center.setX(sum_x/this->noOfPoints);
	center.setY(sum_y/this->noOfPoints);
	center.setZ(sum_z/this->noOfPoints);
};

Point* Face::getCenter(){
	assert(this->center.getId() == this->getId() && "Face::GetCenter() returns error. Calculate the center first.\n");
	return &this->center;
};


void Face::setVertexSign(int vertexNo, int coord, double sign){
	this->vertexSign.setValue(vertexNo, coord, sign);
};

double Face::getVertexSign(int vertexNo, int coord){
	return this->vertexSign.getValue(vertexNo, coord);
};

TensorO2<double>* Face::getVertexSignArray(){
	return &this->vertexSign;
};

void Face::assertVertexSigns(){
	enum coord{x,y,z};
	if (this->getFaceType() == FaceType::Quad){
		assert(this->getVertexSign(0,x)== -1.0 && this->getVertexSign(0,y)== -1.0 );
		assert(this->getVertexSign(1,x)==  1.0 && this->getVertexSign(1,y)== -1.0 );
		assert(this->getVertexSign(2,x)==  1.0 && this->getVertexSign(2,y)==  1.0 );
		assert(this->getVertexSign(3,x)== -1.0 && this->getVertexSign(3,y)==  1.0 );
	};
};


bool Face::findPoint(Point *p){
	bool flag = false;
	assert(this->noOfPoints > 0 && "Face::findPoint() returns error. The face has no points right now.\n");
	for (Index i=0; i<this->getNoOfPoints(); i++){
		if (p == this->getDefiningPoint(i)){
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

double Face::getFlux(int variable, int DOF){
	return this->Flux.getValue(variable, DOF);
};

TensorO2<double>* Face::getFluxArray(){
	return &this->Flux;
};

void Face::setFlux(int variable, int DOF, double Value){
	this->Flux.setValue(variable, DOF, Value);
};

double Face::getPrimitiveVariable(int variable, int DOF){
	return this->primitiveVariableVector.getValue(variable, DOF);
};

TensorO2<double>* Face::getPrimitiveVariableArray(){
	return &this->primitiveVariableVector;
};

void Face::setPrimitiveVariable(int variable, int DOF, double Value){
	this->primitiveVariableVector.setValue(variable, DOF, Value);
};


void Face::setReferenceFace(RefFace* refFace){
	this->refFace = refFace;
	this->refFaceFlag = true;
};

RefFace* Face::getReferenceFace(){
	assert (this->refFaceFlag == true && "Face::getReferenceFace() returns error. Make sure that nDOF number is assigned.\n");
	return this->refFace;
};


void Face::setNDOF(int nDOF){
	this->nDOF = nDOF;
	this->nDOFFlag = true;
};

int Face::getNDOF()const{
	assert (this->nDOFFlag == true && "Face::getNDOF() returns error. Make sure that nDOF number is assigned.\n");
	return this->nDOF;
};

void Face::setNQuad(int nQuad){
	this->nQuad = nQuad;
	this->nQuadFlag = true;
};

int Face::getNQuad()const{
	assert (this->nQuadFlag == true && "Face::getNQuad() returns error. Make sure that nQuad number is assigned.\n");
	return this->nQuad;
};

TensorO1<double>* Face::getJacobian(){
	assert(this->JacobianFlag == true && "Face::getJacobian() returns error. Make sure that the Jacobian is computed.\n");
	return &this->J;
};


TensorO2<double>* Face::getQuadPointsGlobalLocation(){
	return &this->quadPointsGlobalLocation;
};

TensorO2<double>* Face::getDOFPointsGlobalLocation(){
	return &this->DOFPointsGlobalLocation;
};





Matrix<double>* Face::getRotationMatrixParallelToXY(){
	return &this->rotationMatrixParallelToXY;
};

Matrix<double>* Face::getInverseRotationMatrixParallelToXY(){
	return &this->inverseRotationMatrixParallelToXY;
};

void Face::rotateParallelToXY(TensorO1<double> *OldPt, TensorO1<double> *NewPt){
	Math::dot(this->getRotationMatrixParallelToXY(), OldPt, NewPt);
};

void Face::rotateBackParallelToXY(TensorO1<double> *OldPt, TensorO1<double> *NewPt){
	assert(this->getRotationMatrixParallelToXY()->inverseFlag == true);
	Math::dot(this->getInverseRotationMatrixParallelToXY(), OldPt, NewPt);
};

void Face::calculateRotationMatrixParallelToXY(){
	// Calculates Rotation matrix such that it makes normal vector parallel to the z axis (i.e. makes the face parallel to the xy plane)

	// Set size of the rotation matrix;


	// This matrix will be used to rotate the face in order to find the sequence of the points for Jacobian computation
	// Refer Face::calculateJacobianQuadCell() for implementation details.

	// theta: angle of normal vector with the z axis
	// delta: angle of normal vector with the x axis

	assert(this->getRotationMatrixParallelToXY()->getSize0() == 3 and this->getRotationMatrixParallelToXY()->getSize1() == 3);
	assert(Math::approx(this->getNormal()->getL2Norm(), 1.0));

	// first theta is computed:
	double theta = acos(this->normal.getValue(2)); //remember, mag(normal) = 1.0 and nz = n cos(theta)
	//then delta is computed
	double delta = acos(this->normal.getValue(0) / (sin(theta)+FEMTO)); // n*sin(theta)*cos(delta) = nx
	double delta_y = asin(this->normal.getValue(1) / (sin(theta)+FEMTO)); 

	if (this->normal.getValue(1) < 0.0){ //because rotation in xy plane is defined only anticlockwise 
		delta = 2.0*PI - delta; 
	};
	
	//R1: Rotation matrix to rotate normal through delta i.e. get normal in x-z plane
	Matrix<double> R1(3);

	//R2: Rotation matrix to rotate the normal through theta i.e. get normal || to z axis
	Matrix<double> R2(3);

	R1.setValue(0, 0, cos(delta));    R1.setValue(0, 1, sin(delta));   R1.setValue(0, 2, 0.0) ;
	R1.setValue(1, 0, -sin(delta));   R1.setValue(1, 1, cos(delta));   R1.setValue(1, 2, 0.0) ;
	R1.setValue(2, 0, 0.0);           R1.setValue(2, 1, 0.0);          R1.setValue(2, 2, 1.0) ;
	
	R2.setValue(0, 0, cos(theta));    R2.setValue(0, 1, 0.0); 	   R2.setValue(0, 2, -sin(theta)) ;
	R2.setValue(1, 0, 0.0);   	  R2.setValue(1, 1, 1.0);  	   R2.setValue(1, 2, 0.0 );
	R2.setValue(2, 0, sin(theta));    R2.setValue(2, 1, 0.0);  	   R2.setValue(2, 2, cos(theta)) ;


	Math::dot(&R2, &R1, this->getRotationMatrixParallelToXY()); 

	//matrix computation complete


	// confirm whether the matrix is performing okay:
	TensorO1<double> rotNormal(3);
	this->rotateParallelToXY(this->getNormal(), &rotNormal);
	//expected value of the rotNormal is (0,0,1): Let's check that (note: MICRO is pow(10,-6))
	assert(rotNormal.getValue(0) < MICRO    &&     rotNormal.getValue(1) < MICRO      &&     abs(1.0 - rotNormal.getValue(2)) < MICRO     &&      "Please check the rotation matrix computations");
};


void Face::calculateInverseRotationMatrixParallelToXY(){
	assert(this->getInverseRotationMatrixParallelToXY()->getSize0() == 3 and this->getInverseRotationMatrixParallelToXY()->getSize1() == 3);
	this->getRotationMatrixParallelToXY()->getInverse(this->getInverseRotationMatrixParallelToXY());
};





Matrix<double>* Face::getRotationMatrixNormal(){
	return &this->rotationMatrixNormal;
};

Matrix<double>* Face::getInverseRotationMatrixNormal(){
	return &this->inverseRotationMatrixNormal;
};

void Face::rotateNormal(TensorO1<double> *OldPt, TensorO1<double> *NewPt){
	Math::dot(this->getRotationMatrixNormal(), OldPt, NewPt);
};

void Face::rotateBackNormal(TensorO1<double> *OldPt, TensorO1<double> *NewPt){
	assert(this->getRotationMatrixNormal()->inverseFlag == true);
	Math::dot(this->getInverseRotationMatrixNormal(), OldPt, NewPt);
};

void Face::calculateRotationMatrixNormal(){
	// Calculates Rotation matrix such that it projects any given vector in normal-tangent1-tangent2 direction
	assert(this->getRotationMatrixNormal()->getSize0() == 3 and this->getRotationMatrixNormal()->getSize1() == 3);

	for (Index i=0; i<3; i++){
		this->getRotationMatrixNormal()->setValue(0, i, this->getNormal()->getValue(i));
		this->getRotationMatrixNormal()->setValue(1, i, this->getTangent1()->getValue(i));
		this->getRotationMatrixNormal()->setValue(2, i, this->getTangent2()->getValue(i));
	};
	//matrix computation complete
};


void Face::calculateInverseRotationMatrixNormal(){
	assert(this->getInverseRotationMatrixNormal()->getSize0() == 3 and this->getInverseRotationMatrixNormal()->getSize1() == 3);
	this->getRotationMatrixNormal()->getInverse(this->getInverseRotationMatrixNormal());
};






void Face::mapRSTtoXYZQuad(Point *localPt, Point *globalPt){
	// The shape functions are given in: http://textofvideo.nptel.ac.in/105106051/lec40.pdf, page number 8.

	// Check that the vertices have the correct sign 
	this->assertVertexSigns();

	// r,s,t are the coordinates in the local (cell) frame of reference, i.e. defined in the reference cell.
	double r = localPt->getX();
	double s = localPt->getY();
	double t = localPt->getZ(); 	//useless coordinate
	assert(Math::approx(t, 0.0));

	// Enum coord is used for assesing vertexSigns in x or y or z direction. Refer Face::vertexSign and GeometryIO::assignVertexSigns
	// for details.

	enum coord {x,y};

	// XYZ is the local container for global coordinates
	double XYZ[3];
	for (Index i=0; i<3; i++){
		XYZ[i] = 0.0;
	};

	double sign[2];
	for(Index vertex=0; vertex<this->getNoOfPoints(); vertex++){

		sign[0] = this->getVertexSign(vertex, x);
		sign[1] = this->getVertexSign(vertex, y);

		double shape_function_i = 1.0/4.0 * (1.0+sign[0]*r)*(1.0+sign[1]*s);

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
double Face::calculateJacobianQuadFace(double r,double s) {
	
	this->assertVertexSigns();

	double J[3][3];
	// set values to zero (initialization)
	for (Index i=0; i<3; i++){
		for (Index j=0; j<3; j++){
			J[i][j] = 0.0;
		};
	};

	double detJ;
	enum coord{x,y};
	
	for(Index vertex=0; vertex<this->getNoOfPoints(); vertex++){
		double signX = this->getVertexSign(vertex, x);
		double signY = this->getVertexSign(vertex, y);

		J[0][0] += signX * (1.0 + signY * s)  / 4.0 * this->getDefiningPoint(vertex)->getX(); //dx_by_dr
		J[0][1] += signY * (1.0 + signX * r)  / 4.0 * this->getDefiningPoint(vertex)->getX(); //dx_by_ds

		J[1][0] += signX * (1.0 + signY * s)  / 4.0 * this->getDefiningPoint(vertex)->getY(); //dy_by_dr
		J[1][1] += signY * (1.0 + signX * r)  / 4.0 * this->getDefiningPoint(vertex)->getY(); //dy_by_ds

		J[2][0] += signX * (1.0 + signY * s)  / 4.0 * this->getDefiningPoint(vertex)->getZ(); //dz_by_dr
		J[2][1] += signY * (1.0 + signX * r)  / 4.0 * this->getDefiningPoint(vertex)->getZ(); //dz_by_ds

	}
	
	// now the Jacobian matrix is ready
	// det(J) = ||dX/dr x dX/ds||
	double dxbydr = J[0][0]; double dybydr = J[1][0];  double dzbydr = J[2][0]; 
	double dxbyds = J[0][1]; double dybyds = J[1][1];  double dzbyds = J[2][1];

	double term1 = dybydr*dzbyds - dybyds*dzbydr;
	double term2 = -(dxbydr*dzbyds - dxbyds*dzbydr);
	double term3 = dxbydr*dybyds - dxbyds*dybydr;

	detJ = sqrt(term1*term1 + term2*term2 + term3*term3);
	return detJ;
};


int Face::getOwnerQuadPoint(int faceQP){
	return this->mapOwnerQuadPoints.getValue(faceQP);
};

TensorO1<int>* Face::getOwnerQuadPointsArray(){
	return &this->mapOwnerQuadPoints;
};

int Face::getNeighbourQuadPoint(int faceQP){
	return this->mapNeighbourQuadPoints.getValue(faceQP);
};

TensorO1<int>* Face::getNeighbourQuadPointsArray(){
	return &this->mapNeighbourQuadPoints;
};

int Face::getOwnerDOFPoint(int faceDOF){
	return this->mapOwnerDOFPoints.getValue(faceDOF);
};

TensorO1<int>* Face::getOwnerDOFPointsArray(){
	return &this->mapOwnerDOFPoints;
};

int Face::getNeighbourDOFPoint(int faceDOF){
	return this->mapNeighbourDOFPoints.getValue(faceDOF);
};

TensorO1<int>* Face::getNeighbourDOFPointsArray(){
	return &this->mapNeighbourDOFPoints;
};




/*******************************************TESTS**************************************************/


namespace Test{
	void TestFace(){
		cout << "\nRunning test on class Face... " << endl;

		TestFaceInitialization();
		TestDefiningPoints();
		TestNormalAndTangents();
		TestFaceArea();
		TestCenter();

		cout << "Test::TestFace() passed.\n";
	};

	void TestFaceInitialization(){
		int random = getRandom(-10000, 10000);
		Face face(random);
		int grandom = face.getId();
		assert(random == grandom && "Test::TestFaceInitialization() returns error.\n");
	};


	void TestDefiningPoints(){
		int random = getRandom(-10000, 10000);
		Face face(random);

		// create point 1
		double x1, y1, z1; int id1;
		x1 = getRandom(-1000, 1000);
		y1 = getRandom(-1000, 1000);
		z1 = getRandom(-1000, 1000);
		id1 = getRandom(-1000, 1000);
		Point p1(x1, y1, z1, id1);

		// create point 2
		double x2, y2, z2; int id2;
		x2 = getRandom(-2000, 2000);
		y2 = getRandom(-2000, 2000);
		z2 = getRandom(-2000, 2000);
		id2 = getRandom(-2000, 2000);
		Point p2(x2, y2, z2, id2);

		// create point 3
		double x3, y3, z3; int id3;
		x3 = getRandom(-3000, 3000);
		y3 = getRandom(-3000, 3000);
		z3 = getRandom(-3000, 3000);
		id3 = getRandom(-3000, 3000);
		Point p3(x3, y3, z3, id3);

		face.addDefiningPoint(&p1);
		face.addDefiningPoint(&p2);
		face.addDefiningPoint(&p3);

		// test with 3 points 
		unsigned int Np = face.getNoOfPoints();
		assert(Np == 3 && "Test::TestDefiningPoints() in Face.cpp returns error. Face::getNoOfPoints() returns wrong value.\n");

		// create point 4
		double x4, y4, z4; int id4;
		x4 = getRandom(-4000, 4000);
		y4 = getRandom(-4000, 4000);
		z4 = getRandom(-4000, 4000);
		id4 = getRandom(-4000, 4000);
		Point p4(x4, y4, z4, id4);

		// test with 4 points
		face.addDefiningPoint(&p4);
		Np = face.getNoOfPoints();
		assert(Np == 4 && "Test::TestDefiningPoints() in Face.cpp returns error. Face::getNoOfPoints() returns wrong value.\n");


		int idp1 = p1.getId();
		int idp2 = p2.getId();
		int idp3 = p3.getId();
		int idp4 = p4.getId();


		int idf1 = face.getDefiningPoint(0)->getId();
		int idf2 = face.getDefiningPoint(1)->getId();
		int idf3 = face.getDefiningPoint(2)->getId();
		int idf4 = face.getDefiningPoint(3)->getId();

		assert(idp1 == idf1 and idp2 == idf2 and idp3 == idf3 and idp4 == idf4 and "Test::TestDefiningPoints() in Face.cpp returns error. Face::getDefiningPoint() does not return the correct id.\n");

		double xf1, yf1, zf1;
		double xf2, yf2, zf2;
		double xf3, yf3, zf3;
		double xf4, yf4, zf4;

		xf1 = face.getDefiningPoint(0)->getX();
		yf1 = face.getDefiningPoint(0)->getY();
		zf1 = face.getDefiningPoint(0)->getZ();

		xf2 = face.getDefiningPoint(1)->getX();
		yf2 = face.getDefiningPoint(1)->getY();
		zf2 = face.getDefiningPoint(1)->getZ();

		xf3 = face.getDefiningPoint(2)->getX();
		yf3 = face.getDefiningPoint(2)->getY();
		zf3 = face.getDefiningPoint(2)->getZ();

		xf4 = face.getDefiningPoint(3)->getX();
		yf4 = face.getDefiningPoint(3)->getY();
		zf4 = face.getDefiningPoint(3)->getZ();

		assert(xf1 == x1 and yf1 == y1 and zf1 == z1 and
		       xf2 == x2 and yf2 == y2 and zf2 == z2 and
                       xf3 == x3 and yf3 == y3 and zf3 == z3 and
		       "Test::TestDefiningPoints() in Face.cpp returns error. Face::getDefiningPoint() does not return the correct xyz coordinate.\n");


		face.setFaceType();
		FaceType::faceType ft = face.getFaceType();
		assert (ft == FaceType::Quad && "Test::TestDefiningPoints() in Face.cpp returns error. FaceType based on the defining points not set correctly.\n");


	};

	void TestNormalAndTangents(){
		int random = getRandom(-10000, 10000);

		// create point 1
		double x1, y1, z1; int id1;
		x1 = getRandom(-1000.0, 1000.0);
		y1 = getRandom(-1000.0, 1000.0);
		z1 = getRandom(-1000.0, 1000.0);
		id1 = getRandom(-1000.0, 1000.0);
		Point p1(x1, y1, z1, id1);

		// create point 2
		double x2, y2, z2; int id2;
		x2 = getRandom(-2000.0, 2000.0);
		y2 = getRandom(-2000.0, 2000.0);
		z2 = getRandom(-2000.0, 2000.0);
		id2 = getRandom(-2000.0, 2000.0);
		Point p2(x2, y2, z2, id2);

		// create point 3
		double x3, y3, z3; int id3;
		x3 = getRandom(-3000.0, 3000.0);
		y3 = getRandom(-3000.0, 3000.0);
		z3 = getRandom(-3000.0, 3000.0);
		id3 = getRandom(-3000.0, 3000.0);
		Point p3(x3, y3, z3, id3);

		// create point 4
		double x4, y4, z4; int id4;
		x4 = getRandom(-4000.0, 4000.0);
		y4 = getRandom(-4000.0, 4000.0);
		z4 = getRandom(-4000.0, 4000.0);
		id4 = getRandom(-4000.0, 4000.0);
		Point p4(x4, y4, z4, id4);

		// Testing for a tri face
		Face faceTri(random);
		// Testing for a quad face
		Face faceQuad(random);

		faceTri.addDefiningPoint(&p1);
		faceTri.addDefiningPoint(&p2);
		faceTri.addDefiningPoint(&p3);

		faceQuad.addDefiningPoint(&p1);
		faceQuad.addDefiningPoint(&p2);
		faceQuad.addDefiningPoint(&p3);
		faceQuad.addDefiningPoint(&p4);


		faceTri.setFaceType();
		faceTri.calculateNormal();

		double normalLength = faceTri.getNormal()->getL2Norm();
		assert (abs(normalLength - 1.0) < PICO && "Test::TestNormalAndTangents() in Face.cpp returns error. Computed normal does not have a unit length.\n");

		double tan1Length = faceTri.getTangent1()->getL2Norm();
		assert (abs(tan1Length - 1.0) < PICO && "Test::TestNormalAndTangents() in Face.cpp returns error. Computed tangent1 does not have a unit length.\n");

		double tan2Length = faceTri.getTangent2()->getL2Norm();
		assert (abs(tan2Length - 1.0) < PICO && "Test::TestNormalAndTangents() in Face.cpp returns error. Computed tangent2 does not have a unit length.\n");

		// make sure that the Tangents are computed correctly
		assert(abs(Math::dot(faceTri.getNormal(),faceTri.getTangent1())) < NANO && "Face::calculateNormal() returns error. Check the Tangent vector computation in Face.cpp");
		assert(abs(Math::dot(faceTri.getTangent1(),faceTri.getTangent2())) < NANO && "Face::calculateNormal() returns error. Check the Tangent vector computation in Face.cpp");
		assert(abs(Math::dot(faceTri.getNormal(),faceTri.getTangent2())) < NANO && "Face::calculateNormal() returns error. Check the Tangent vector computation in Face.cpp");




		faceQuad.setFaceType();
		faceQuad.calculateNormal();

		normalLength = faceQuad.getNormal()->getL2Norm();
		assert (abs(normalLength - 1.0) < PICO && "Test::TestNormalAndTangents() in Face.cpp returns error. Computed normal does not have a unit length.\n");

		tan1Length = faceQuad.getTangent1()->getL2Norm();
		assert (abs(tan1Length - 1.0) < PICO && "Test::TestNormalAndTangents() in Face.cpp returns error. Computed tangent1 does not have a unit length.\n");

		tan2Length = faceQuad.getTangent2()->getL2Norm();
		assert (abs(tan2Length - 1.0) < PICO && "Test::TestNormalAndTangents() in Face.cpp returns error. Computed tangent2 does not have a unit length.\n");

		// make sure that the Tangents are computed correctly
		assert(abs(Math::dot(faceQuad.getNormal(),faceQuad.getTangent1())) < NANO && "Face::calculateNormal() returns error. Check the Tangent vector computation in Face.cpp");
		assert(abs(Math::dot(faceQuad.getTangent1(),faceQuad.getTangent2())) < NANO && "Face::calculateNormal() returns error. Check the Tangent vector computation in Face.cpp");
		assert(abs(Math::dot(faceQuad.getNormal(),faceQuad.getTangent2())) < NANO && "Face::calculateNormal() returns error. Check the Tangent vector computation in Face.cpp");


	};

	void TestFaceArea(){
		double x1, y1, z1; int id1;
		double x2, y2, z2; int id2;
		double x3, y3, z3; int id3;
		double x4, y4, z4; int id4;


		double areaTriA ; 
                double areaQuadA ;
                                  
                double ToleranceQ;
                double ToleranceT;

		for (Index i=0; i<1000; i++){
			// Testing for regular points (i.e. well arranged)  

			// 1. Points arranged counterclockwise
			int random = getRandom(-10000, 10000);

			// create point 1 in the first quadrant
			x1 = getRandom(0.0, 50.0);
			y1 = getRandom(50.0, 100.0);
			z1 = random;
			id1 = getRandom(0, 1000);
			Point p1(x1, y1, z1, id1);

			// create point 2 in the second quadrant
			x2 = getRandom(-100.0, -50.0);
			y2 = getRandom(0.0, 50.0);
			z2 = random;
			id2 = getRandom(0, 1000);
			Point p2(x2, y2, z2, id2);

			// create point 3 in the third quadrant
			x3 = getRandom(-50.0, 0.0);
			y3 = getRandom(-100.0, -50.0);
			z3 = random;
			id3 = getRandom(0, 1000);
			Point p3(x3, y3, z3, id3);

			// create point 4 in the fourth quadrant
			x4 = getRandom(50.0, 100.0);
			y4 = getRandom(-50.0, 0.0);
			z4 = random;
			id4 = getRandom(0, 1000);
			Point p4(x4, y4, z4, id4);

			// Testing for a tri face
			Face faceTri(random);
			// Testing for a quad face
			Face faceQuad(random);

			faceTri.addDefiningPoint(&p1);
			faceTri.addDefiningPoint(&p2);
			faceTri.addDefiningPoint(&p3);

			faceQuad.addDefiningPoint(&p1);
			faceQuad.addDefiningPoint(&p2);
			faceQuad.addDefiningPoint(&p3);
			faceQuad.addDefiningPoint(&p4);


			faceTri.setFaceType();
			faceTri.calculateNormal();
			faceTri.calculateArea();

			faceQuad.setFaceType();
			faceQuad.calculateNormal();
			faceQuad.calculateArea();

			areaTriA = Math::calculateAreaOfTriangle(p1, p2, p3);
			areaQuadA = Math::calculateAreaOfTriangle(p1, p2, p3) + Math::calculateAreaOfTriangle(p4, p1, p3);

			ToleranceQ = areaQuadA * 0.001;
			ToleranceT = areaTriA * 0.001;

			assert(abs(areaTriA - faceTri.getArea()) < ToleranceT && "Test::TestArea() in Face.cpp returns error. Make sure that area for a triangle is computed correctly.\n");
			assert(abs(areaQuadA - faceQuad.getArea()) < ToleranceQ && "Test::TestArea() in Face.cpp returns error. Make sure that area for a quadrilateral is computed correctly.\n");

		};

		for (Index i=0; i<1000; i++){
			// Testing for regular points (i.e. well arranged)  

			// 2. Points arranged clockwise
			int random = getRandom(-10000, 10000);

			// create point 1 in the second quadrant
			x1 = getRandom(-100.0, -50.0);
			y1 = getRandom(0.0, 50.0);
			z1 = random;
			id1 = getRandom(0, 1000);
			Point p1(x1, y1, z1, id1);

			// create point 2 in the first quadrant
			x2 = getRandom(0.0, 50.0);
			y2 = getRandom(50.0, 100.0);
			z2 = random;
			id2 = getRandom(0, 1000);
			Point p2(x2, y2, z2, id2);

			// create point 3 in the fourth quadrant
			x3 = getRandom(50.0, 100.0);
			y3 = getRandom(-50.0, 0.0);
			z3 = random;
			id3 = getRandom(0, 1000);
			Point p3(x3, y3, z3, id3);

			// create point 4 in the third quadrant
			x4 = getRandom(-50.0, 0.0);
			y4 = getRandom(-100.0, -50.0);
			z4 = random;
			id4 = getRandom(0, 1000);
			Point p4(x4, y4, z4, id4);


			// Testing for a tri face
			Face faceTri1(random);
			// Testing for a quad face
			Face faceQuad1(random);

			faceTri1.addDefiningPoint(&p1);
			faceTri1.addDefiningPoint(&p2);
			faceTri1.addDefiningPoint(&p3);

			faceQuad1.addDefiningPoint(&p1);
			faceQuad1.addDefiningPoint(&p2);
			faceQuad1.addDefiningPoint(&p3);
			faceQuad1.addDefiningPoint(&p4);


			faceTri1.setFaceType();
			faceTri1.calculateNormal();
			faceTri1.calculateArea();

			faceQuad1.setFaceType();
			faceQuad1.calculateNormal();
			faceQuad1.calculateArea();

			areaTriA = Math::calculateAreaOfTriangle(p1, p2, p3);
			areaQuadA = Math::calculateAreaOfTriangle(p1, p2, p3) + Math::calculateAreaOfTriangle(p4, p1, p3);

			ToleranceQ = areaQuadA * 0.001;
			ToleranceT = areaTriA * 0.001;

			assert(abs(areaTriA - faceTri1.getArea()) < ToleranceT && "Test::TestArea() in Face.cpp returns error. Make sure that area for a triangle is computed correctly.\n");
			assert(abs(areaQuadA - faceQuad1.getArea()) < ToleranceQ && "Test::TestArea() in Face.cpp returns error. Make sure that area for a quadrilateral is computed correctly.\n");
		};

		for (Index i=0; i<1000; i++){
			//3. Testing for irregular points  

			int random = getRandom(-10000, 10000);

			// create point 1 in the second quadrant
			x1 = getRandom(-100.0, -50.0);
			y1 = getRandom(0.0, 50.0);
			z1 = random;
			id1 = getRandom(0, 1000);
			Point p1(x1, y1, z1, id1);

			// create point 2 in the fourth quadrant
			x2 = getRandom(50.0, 100.0);
			y2 = getRandom(-50.0, 0.0);
			z2 = random;
			id2 = getRandom(0, 1000);
			Point p2(x2, y2, z2, id2);

			// create point 3 in the first quadrant
			x3 = getRandom(0.0, 50.0);
			y3 = getRandom(50.0, 100.0);
			z3 = random;
			id3 = getRandom(0, 1000);
			Point p3(x3, y3, z3, id3);

			// create point 4 in the third quadrant
			x4 = getRandom(-50.0, 0.0);
			y4 = getRandom(-100.0, -50.0);
			z4 = random;
			id4 = getRandom(0, 1000);
			Point p4(x4, y4, z4, id4);


			// Testing for a tri face
			Face faceTri1(random);
			// Testing for a quad face
			Face faceQuad1(random);

			faceTri1.addDefiningPoint(&p1);
			faceTri1.addDefiningPoint(&p2);
			faceTri1.addDefiningPoint(&p3);

			faceQuad1.addDefiningPoint(&p1);
			faceQuad1.addDefiningPoint(&p2);
			faceQuad1.addDefiningPoint(&p3);
			faceQuad1.addDefiningPoint(&p4);


			faceTri1.setFaceType();
			faceTri1.calculateNormal();
			faceTri1.calculateArea();

			faceQuad1.setFaceType();
			faceQuad1.calculateNormal();
			faceQuad1.calculateArea();

			areaTriA = Math::calculateAreaOfTriangle(p1, p2, p3);
			areaQuadA = Math::calculateAreaOfTriangle(p1, p2, p3) + Math::calculateAreaOfTriangle(p4, p1, p2);

			ToleranceQ = areaQuadA * 0.001;
			ToleranceT = areaTriA * 0.001;

			assert(abs(areaTriA - faceTri1.getArea()) < ToleranceT && "Test::TestArea() in Face.cpp returns error. Make sure that area for a triangle is computed correctly.\n");
			assert(abs(areaQuadA - faceQuad1.getArea()) < ToleranceQ && "Test::TestArea() in Face.cpp returns error. Make sure that area for a quadrilateral is computed correctly.\n");
		};


		for (Index i=0; i<1000; i++){
			//4. Testing for rectangles (area compared analytically)  

			int random = getRandom(-10000, 10000);

			// create point 1 in the second quadrant
			x1 = getRandom(-100.0, -50.0);
			y1 = getRandom(0.0, 50.0);
			z1 = random;
			id1 = getRandom(0, 1000);
			Point p1(x1, y1, z1, id1);

			// create point 2 in the fourth quadrant
			x2 = getRandom(50.0, 100.0);
			y2 = p1.getY();
			z2 = random;
			id2 = getRandom(0, 1000);
			Point p2(x2, y2, z2, id2);

			double length = abs(x2 - x1);

			// create point 3 in the first quadrant
			x3 = x2;
			y3 = getRandom(-100.0, 100.0);
			z3 = random;
			id3 = getRandom(0, 1000);
			Point p3(x3, y3, z3, id3);

			double height = abs(y3 - y2);

			// create point 4 in the third quadrant
			x4 = x1;
			y4 = y3;
			z4 = random;
			id4 = getRandom(0, 1000);
			Point p4(x4, y4, z4, id4);


			// Testing for a quad face
			Face faceQuad1(random);

			faceQuad1.addDefiningPoint(&p1);
			faceQuad1.addDefiningPoint(&p2);
			faceQuad1.addDefiningPoint(&p3);
			faceQuad1.addDefiningPoint(&p4);


			faceQuad1.setFaceType();
			faceQuad1.calculateNormal();
			faceQuad1.calculateArea();

			areaQuadA = length * height;

			ToleranceQ = NANO;

			assert(abs(areaQuadA - faceQuad1.getArea()) < ToleranceQ && "Test::TestArea() in Face.cpp returns error. Make sure that area for a quadrilateral is computed correctly.\n");
		};

	};

	
	void TestCenter(){
		double x1, y1, z1; int id1;
		double x2, y2, z2; int id2;
		double x3, y3, z3; int id3;
		double x4, y4, z4; int id4;
		double xC, yC, zC;

		for (Index i=0; i<1000; i++){
			//1. Testing for rectangles  

			int random = getRandom(-10000, 10000);

			// create point 1 in the second quadrant
			x1 = getRandom(-100.0, -50.0);
			y1 = getRandom(0.0, 50.0);
			z1 = random;
			id1 = getRandom(0, 1000);
			Point p1(x1, y1, z1, id1);

			// create point 2 in the fourth quadrant
			x2 = getRandom(50.0, 100.0);
			y2 = p1.getY();
			z2 = random;
			id2 = getRandom(0, 1000);
			Point p2(x2, y2, z2, id2);

			double length = abs(x2 - x1);

			// create point 3 in the first quadrant
			x3 = x2;
			y3 = getRandom(-100.0, 100.0);
			z3 = random;
			id3 = getRandom(0, 1000);
			Point p3(x3, y3, z3, id3);

			double height = abs(y3 - y2);

			// create point 4 in the third quadrant
			x4 = x1;
			y4 = y3;
			z4 = random;
			id4 = getRandom(0, 1000);
			Point p4(x4, y4, z4, id4);


			// Testing for a quad face
			Face faceQuad1(random);

			faceQuad1.addDefiningPoint(&p1);
			faceQuad1.addDefiningPoint(&p2);
			faceQuad1.addDefiningPoint(&p3);
			faceQuad1.addDefiningPoint(&p4);


			faceQuad1.setFaceType();
			faceQuad1.calculateCenter();

			xC = (x1 + x2 + x3 + x4)/4.0;
			yC = (y1 + y2 + y3 + y4)/4.0;
			zC = (z1 + z2 + z3 + z4)/4.0;

			double ToleranceQ = NANO;

			assert(faceQuad1.getCenter()->getX() == xC && "Test::TestCenter() in Face.cpp returns error. Make sure that center for a quadrilateral is computed correctly.\n");
			assert(faceQuad1.getCenter()->getY() == yC && "Test::TestCenter() in Face.cpp returns error. Make sure that center for a quadrilateral is computed correctly.\n");
			assert(faceQuad1.getCenter()->getZ() == zC && "Test::TestCenter() in Face.cpp returns error. Make sure that center for a quadrilateral is computed correctly.\n");
		};
	};
};

