#include "RefCell.h"


//_________________________________________________Base class__________________________________________________//
// Reference cell

//constructor
RefCell :: RefCell(){
};

//destructor
RefCell :: ~RefCell(){
};

int RefCell :: getNoOfPoints()const{
	return this->noOfPoints;
};

void RefCell :: setNoOfPoints(int pts){
	this->noOfPoints = pts;
};

int RefCell :: getNoOfFaces()const{
	return this->noOfFaces;
};

void RefCell :: setNoOfFaces(int fcs){
	this->noOfFaces = fcs;
};

int RefCell :: getNoOfDOFPoints()const{
	return this->nDOF;
};

void RefCell :: setNoOfDOFPoints(int dofs){
	this->nDOF = dofs;
};

void RefCell::setIntegrationType(IntFlag::intflag intFlag){
	assert(intFlag == IntFlag::inexact or intFlag == IntFlag::exact && "Set integration type to inexact (intFlag 0) or exact (intFlag 1)");
	this->intFlag = intFlag;
};


TensorO1<double>* RefCell::getQuadWeights(){
	return &this->w_q;
};


void RefCell::init(int order, IntFlag::intflag intFlag, CellType::cellType celltype) {
	// order: polynomial order. 
	// intFlag: integration type. 0:inexact, 1:exact
	assert(intFlag == 0 or intFlag == 1 && "Set integration type to inexact (intFlag 0) or exact (intFlag 1)");

	// Set the cell type
	this->celltype = celltype; 
	// Set the order
	this->order = order;
	// Set the integration flag (exact vs inexact)
	this->intFlag = intFlag;

	assert (celltype == CellType::Hex or celltype == CellType::Tet or celltype == CellType::Prism or celltype == CellType::Pyramid);

	if (celltype == CellType::Hex){
		this->noOfPoints = 8;
		this->noOfFaces = 6;
		this->nDOF = pow((this->order + 1),3); //coded for a 3d goemetry only.

		int K = pow((this->order + 1 + this->intFlag),3);

		//setting up the cell Jacobian at quadrature points. For a refCell, all det(J) = 1.0
		J.setSize(K);

		for (int i=0; i<K; i++){
			J.setValue(i, 1.0);
		};
	}

	else if (celltype == CellType::Tet){
	}
	else if (celltype == CellType::Prism){
	}
	else if (celltype == CellType::Pyramid){
	};

	this->generateQuadratureWeights();
	this->generateLagMatrix();
	this->generateGradLagMatrix();
}


CellType::cellType RefCell::getCellType(){
	return this->celltype;
};

void RefCell::generateQuadratureWeights(){
	// Fills array w_q
	FunctionalSpace F(this->order, this->intFlag);
	if (celltype == CellType::Hex){
		int K = this->order + 1 + this->intFlag;
		int nQuad = K*K*K; //total number of quadrature points.
		this->nQuad = nQuad;

		w_q.setSize(nQuad);
		TensorO1<double> x_dummy(nQuad); 
		TensorO1<double> y_dummy(nQuad); 
		TensorO1<double> z_dummy(nQuad);

		F.LGLRootsAndWeights3D(K, K, K, &x_dummy, &y_dummy, &z_dummy, &this->w_q);
			
	}
	else if (celltype == CellType::Tet){
	}
	else if (celltype == CellType::Prism){
	}
	else if (celltype == CellType::Pyramid){
	}
};



void RefCell::generateLagMatrix(){
	// LagMatrix is a matrix of the Lagrange polynomials computed at the quadrature points.
	// It has dim: nQuad x nDOF

	FunctionalSpace F(this->order, this->intFlag);
	if (celltype == CellType::Hex){
		F.generateLagMatrix3DTensor(&this->LagMatrix);
	}
	else if (celltype == CellType::Tet){
		//insert function for generating Lagrange matrix for tet cells
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
	else if (celltype == CellType::Prism){
		//insert function for generating Lagrange matrix for prism
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
	else if (celltype == CellType::Pyramid){
		//insert function for generating Lagrange matrix for pyramid
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
};



void RefCell::generateGradLagMatrix(){
	// LagDrstMatrix is a matrix of the gradient of Lagrange polynomials in r, s, t direction computed at the quadrature points.
	// It has dim: nQuad x nDOF 

	FunctionalSpace F(this->order, this->intFlag);
	if (celltype == CellType::Hex){
		F.generateGradLagMatrix3DTensor(&this->LagDrMatrix, &this->LagDsMatrix, &this->LagDtMatrix);
	}
	else if (celltype == CellType::Tet){
		//insert function for generating Lagrange matrix for tet
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
	else if (celltype == CellType::Prism){
		//insert function for generating Lagrange matrix for prism
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
	else if (celltype == CellType::Pyramid){
		//insert function for generating Lagrange matrix for pyramid
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
};


// Get pointer to the Lagrange matrix
TensorO2<double>* RefCell::getLagMatrix(){
	return &this->LagMatrix;
};


// Get pointer to the Lagrange derivative matrix in x direction (r in reference frame)
TensorO2<double>* RefCell::getLagDrMatrix(){
	return &this->LagDrMatrix;
};

// Get pointer to the Lagrange derivative matrix in y direction (s in reference frame)
TensorO2<double>* RefCell::getLagDsMatrix(){
	return &this->LagDsMatrix;
};

// Get pointer to the Lagrange derivative matrix in z direction (t in reference frame)
TensorO2<double>* RefCell::getLagDtMatrix(){
	return &this->LagDtMatrix;
};



// Other printing functions

void RefCell::printLagMatrix(){
	cout <<"\nLagMat: \n";
	for (int i=0; i<nQuad; i++){
		for (int j=0; j<nDOF; j++){
			cout << LagMatrix.getValue(i,j) << "\t";
		};
		cout << endl;
	};
};


void RefCell::printGradLagMatrix(){
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
	cout <<"\nGradLagMat t: \n";
	for (int i=0; i<nQuad; i++){
		for (int j=0; j<nDOF; j++){
			cout << LagDtMatrix.getValue(i,j) << "\t";
		};
		cout << endl;
	};
};

void RefCell::print(){
	string celltypes[4] = {"Hex", "Tet", "Prism", "Pyramid"}; 

	cout << "\nThe cell type is: " << celltypes[this->celltype] << endl;
	cout << "\nTotal number of defining points are: " << this->noOfPoints << endl;
	cout << "\nTotal number of defining faces are: " << this->noOfFaces << endl;
	cout << "\nOrder of polynomial reconstruction: " << this->order << endl;
	cout << endl;
	cout << "Quadrature weights:  [ ";
	for (int i=0; i<nQuad; i++){
		cout << w_q.getValue(i) << "  " ;
	};
	cout << " ]\n";
};



//*******************************************************************************************************************//
// Tests:

namespace Test{
	void TestReferenceCells(){

		cout << "\nRunning tests on the Reference Cells.\n";

		for (Index i=1; i<3; i++){
			RefCell hexCell;
			hexCell.init(i, IntFlag::inexact, CellType::Hex);
			hexCell.generateQuadratureWeights();

			assert(hexCell.getNoOfPoints() == 8);
			assert(hexCell.getNoOfFaces() == 6);
			assert(hexCell.getNoOfDOFPoints() == pow((i+1),3));
			assert(Math::approx(hexCell.getQuadWeights()->getL1Norm(), 8.0));

		};

		cout << "Test::TestReferenceCells() passed.\n";
	};
};

