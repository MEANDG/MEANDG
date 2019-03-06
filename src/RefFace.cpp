#include "RefFace.h"


// Constructor
RefFace :: RefFace(){
};

//destructor
RefFace::~RefFace(){
};

int RefFace :: getNoOfPoints()const{
	return this->noOfPoints;
};

void RefFace :: setNoOfPoints(int pts){
	this->noOfPoints = pts;
};


int RefFace :: getNoOfDOFPoints()const{
	return this->nDOF;
};

void RefFace :: setNoOfDOFPoints(int dofs){
	this->nDOF = dofs;
};

void RefFace::setIntegrationType(IntFlag::intflag intFlag){
	assert(intFlag == IntFlag::inexact or intFlag == IntFlag::exact && "Set integration type to inexact (intFlag 0) or exact (intFlag 1)");
	this->intFlag = intFlag;
};


TensorO1<double>* RefFace::getQuadWeights(){
	return &this->w_q;
};


void RefFace::init(int order, IntFlag::intflag intFlag, FaceType::faceType facetype){
	// order: polynomial order. 
	// intFlag: integration type. 0:inexact, 1:exact

	assert(intFlag == 0 or intFlag == 1 && "Set integration type to inexact (intFlag 0) or exact (intFlag 1)");

	// Set the face type
	this->facetype = facetype; 
	// Set the order
	this->order = order;
	// Set the integration flag (exact vs inexact)
	this->intFlag = intFlag;

	assert (facetype == FaceType::Tri or facetype == FaceType::Quad);

	if (facetype == FaceType::Quad){
		this->noOfPoints = 4;
		this->nDOF = pow((this->order + 1 ),2);

		int K = pow((this->order + 1 + this->intFlag),2);

		J. setSize(K);
		w_q.setSize(K);

		for (int i=0; i<K ; i++){
			J.setValue(i, 1.0);
			w_q.setValue(i, 0.0); //initialization. Will be filled later.
		};
	}
	else if (facetype == FaceType::Tri){
		// set up the tri face
	};

	// Quadrature weights computed here
	this->generateQuadratureWeights();
	// LagMat generated here. LagMat has nQuad rows and nDOF columns. It stores nDOF Lagrange polynomials at each quadrature point.
	this->generateLagMatrix();
	// Gradient of LagMat (in r and s axis)
	this->generateGradLagMatrix();

};

FaceType::faceType RefFace::getFaceType(){
	return this->facetype;
};

void RefFace::generateQuadratureWeights(){
	// Fills array w_q
	FunctionalSpace F(this->order, this->intFlag);
	if (facetype == FaceType::Quad){
		int K = this->order + 1 + this->intFlag;
		int nQuad = K*K; //total number of quadrature points.
		this->nQuad = nQuad;

		w_q.setSize(nQuad);

		TensorO1<double> x_dummy(nQuad); 
		TensorO1<double> y_dummy(nQuad); 

		F.LGLRootsAndWeights2D(K, K, &x_dummy, &y_dummy, &this->w_q);
		// w_q gets filled here.
	}
	else if (facetype == FaceType::Tri){
	};
};

void RefFace::generateLagMatrix(){
	// LagMatrix is a matrix of the Lagrange polynomials computed at the quadrature points.
	// It has dim: nQuad x nDOF

	FunctionalSpace F(this->order, this->intFlag);

	if (facetype == FaceType::Quad){
		F.generateLagMatrix2DTensor(&LagMatrix); //assumed LGL roots and weights for now
	}
	else if (facetype == FaceType::Tri){
		// Insert function for generating Lagrange matrix for Tri Faces
	};
};


void RefFace::generateGradLagMatrix(){
	// LagDrstMatrix is a matrix of the gradient of Lagrange polynomials in r, s, t direction computed at the quadrature points.
	// It has dim: nQuad x nDOF 

	FunctionalSpace F(this->order, this->intFlag);

	if (facetype == FaceType::Quad){
		F.generateGradLagMatrix2DTensor(&LagDrMatrix, &LagDsMatrix); //assumed LGL roots and weights for now
	}
	else if (facetype == FaceType::Tri){
		// Insert functions for generating Lagrange derivative matrix for Tri face
	}
};

// Get pointer to the Lagrange matrix
TensorO2<double>* RefFace::getLagMatrix(){
	return &this->LagMatrix;
};


// Get pointer to the Lagrange derivative matrix in x direction (r in reference frame)
TensorO2<double>* RefFace::getLagDrMatrix(){
	return &this->LagDrMatrix;
};

// Get pointer to the Lagrange derivative matrix in y direction (s in reference frame)
TensorO2<double>* RefFace::getLagDsMatrix(){
	return &this->LagDsMatrix;
};

// Other printing functions

void RefFace::printLagMatrix(){
	cout <<"\nLagMat: \n";
	for (int i=0; i<nQuad; i++){
		for (int j=0; j<nDOF; j++){
			cout << LagMatrix.getValue(i,j) << "\t";
		};
		cout << endl;
	};
};


void RefFace::printGradLagMatrix(){
	cout <<"\nGradLagMat r: \n";
	for (int i=0; i<nQuad; i++){
		for (int j=0; j<nDOF; j++){
			cout << LagDrMatrix.getValue(i,j) << "\t";
		};
		cout << endl;
	};
	cout <<"\nGradLagMat s: \n";
	for (int i=0; i<nQuad; i++){
		for (int j=0; j<nDOF; j++){
			cout << LagDsMatrix.getValue(i,j) << "\t";
		};
		cout << endl;
	};
};


void RefFace::print(){
	string facetypes[2] = {"quadrilateral","triangle"}; 

	cout << "\nThe face type is: " << facetypes[this->facetype] << endl;
	cout << "\nTotal number of defining points are: " << this->noOfPoints << endl;
	cout << "\nOrder of polynomial reconstruction: " << this->order << endl;
};

