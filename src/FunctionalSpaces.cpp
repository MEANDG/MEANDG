#include "FunctionalSpaces.h"


FunctionalSpace::FunctionalSpace(int order, IntFlag::intflag intFlag){
	assert(order >= 1);
	this->order = order;
	this->intFlag = intFlag;
};

FunctionalSpace::~FunctionalSpace(){
};




unsigned int FunctionalSpace::getOrder()const{
	return this->order;
};

void FunctionalSpace::setOrder(unsigned int order){
	assert(order >= 1);
	this->order = order;
};

void FunctionalSpace::setIntFlag(IntFlag::intflag intFlag){
	this->intFlag = intFlag;
};

IntFlag::intflag FunctionalSpace::getIntFlag(){
	return this->intFlag;
};

//************************************************************************************************//
// Nodal Space
//************************************************************************************************//

// Pass by pointer version
void FunctionalSpace::LagrangePolynomial1D(double r, TensorO1<double> *x, TensorO1<double> *phi){
 	// The function fills the array of the N+1 Lagrange polynomials phi.
	// r: location at which interpolation is to be performed.      
	// x: An array of points where sampling of the data is done. Total number of elements of x = order+1
	//    These points can be equi-spaced or roots of orthonormal polynomials. 
	// phi: array of the Lagrange polynomials phi[i] gives ith Lagrange polynomial

	int Np = this->getOrder() + 1; // for 1D interpolation N+1 polynomials required for Nth order interpolation 

	for (Index i=0; i< Np; i++){
		// computation of ith Lagrange polynomial
		double prod = 1.0;
		for (Index j=0; j<Np ; j++){
			if (j != i){
				prod = prod * (r - x->getValue(j))/(x->getValue(i) - x->getValue(j));		
			};
		};
		phi->setValue(i, prod); //ith Lagrange polynomial computed
	};
};


// Pass by reference version
void FunctionalSpace::LagrangePolynomial1D(double r, TensorO1<double> &x, TensorO1<double> &phi){
 	// The function fills the array of the N+1 Lagrange polynomials phi.
	// r: location at which interpolation is to be performed.      
	// x: An array of points where sampling of the data is done. Total number of elements of x = order+1
	//    These points can be equi-spaced or roots of orthonormal polynomials. 
	// phi: array of the Lagrange polynomials phi[i] gives ith Lagrange polynomial

	int Np = this->getOrder() + 1; // for 1D interpolation N+1 polynomials required for Nth order interpolation 

	for (Index i=0; i< Np; i++){
		// computation of ith Lagrange polynomial
		double prod = 1.0;
		for (Index j=0; j<Np ; j++){
			if (j != i){
				prod = prod * (r - x.getValue(j))/(x.getValue(i) - x.getValue(j));		
			};
		};
		phi.setValue(i, prod); //ith Lagrange polynomial computed
	};
};


// Pass by reference version
void FunctionalSpace::LagrangePolynomialDerivative1D(double r, TensorO1<double> &x, TensorO1<double> &phi){
 	// The function fills the array of derivatives of the N+1 Lagrange polynomials phi. i.e. phi'
	// r: location at which interpolation is to be performed.      
	// x: An array of points where sampling of the data is done. Total number of elements of x = N+1
	//    These points can be equi-spaced or roots of orthonormal polynomials. 
	// phi: array of the Lagrange polynomials phi[i] gives ith Lagrange polynomial

	int Np = this->getOrder() + 1; // for 1D interpolation N+1 polynomials required for Nth order interpolation 

	for (Index i=0; i< Np; i++){
		// computation of derivative of the ith Lagrange polynomial
		double sum = 0;
		for (Index k=0; k<Np; k++){
			if (k!= i){
				double prod = 1.0/(x.getValue(i) - x.getValue(k));
				for (Index j=0; j<Np ; j++){
					if (j != i and j!=k){
						prod = prod * (r - x.getValue(j))/(x.getValue(i) - x.getValue(j));		
					};
				};
				sum += prod; //ith Lagrange polynomial computed
			};
		};
		phi.setValue(i, sum);
	};
};


// Pass by pointer version
void FunctionalSpace::LagrangePolynomialDerivative1D(double r, TensorO1<double> *x, TensorO1<double> *phi){
 	// The function fills the array of derivatives of the N+1 Lagrange polynomials phi. i.e. phi'
	// r: location at which interpolation is to be performed.      
	// x: An array of points where sampling of the data is done. Total number of elements of x = N+1
	//    These points can be equi-spaced or roots of orthonormal polynomials. 
	// phi: array of the Lagrange polynomials phi[i] gives ith Lagrange polynomial

	int Np = this->getOrder() + 1; // for 1D interpolation N+1 polynomials required for Nth order interpolation 

	for (Index i=0; i< Np; i++){
		// computation of derivative of the ith Lagrange polynomial
		double sum = 0;
		for (Index k=0; k<Np; k++){
			if (k!= i){
				double prod = 1.0/(x->getValue(i) - x->getValue(k));
				for (Index j=0; j<Np ; j++){
					if (j != i and j!=k){
						prod = prod * (r - x->getValue(j))/(x->getValue(i) - x->getValue(j));		
					};
				};
				sum += prod; //ith Lagrange polynomial computed
			};
		};
		phi->setValue(i, sum);
	};
};



void FunctionalSpace::LagrangePolynomial2DTensor(double x, double y, TensorO1<double> *xi, TensorO1<double> *phi){
 	// The function fills the array of the N+1 Lagrange polynomials phi.
	// x,y: Location at which interpolation is to be performed.      
	// xi: An array of points where sampling of the data is done. Total number of elements of xi = order+1
	//    These points can be equi-spaced or roots of orthonormal polynomials. 
	// phi: Array of the Lagrange polynomials phi[i] gives ith Lagrange polynomial. 

	int Np = xi->getSize();
	TensorO1<double> phi1Dx(Np); // Npx Lagrange polynomials at location x
	TensorO1<double> phi1Dy(Np); // Npy Lagrange polynomials at location y

	this->LagrangePolynomial1D(x,xi,&phi1Dx); // Lagrange polynomials are computed for location x
	this->LagrangePolynomial1D(y,xi,&phi1Dy); // Lagrange polynomials are computed for location y

	int counter = 0;
	for (Index i=0;i<Np;i++){
		for (Index j=0;j<Np;j++){
			// The Lagrange polynomials are obtained by taking a tensor product of 1D Lagrange polynomials
			phi->setValue(counter, phi1Dx.getValue(i)*phi1Dy.getValue(j));
			counter++;
		}
	}
};


void FunctionalSpace::LagrangePolynomialDerivative2DTensor(double x, double y, TensorO1<double> *xi, TensorO1<double> *phi_x, TensorO1<double> *phi_y){
 	//The function fills derivatives of Lagrange Tensor polynomial in x and in y direction
	// (x,y): location at which differentiation is to be performed.      
	// (xi): An array of points where sampling of the data is done. Total number of elements of xi = N+1
	//    These points can be equi-spaced or roots of orthonormal polynomials. 
	// phi_x: array of the partial derivatives of Lagrange polynomials in x direction. ph_x[i] gives ith component of d(phi)/dx
	// phi_y: array of the partial derivatives of Lagrange polynomials in y direction. ph_y[i] gives ith component of d(phi)/dy
	
	int Np = xi->getSize();
	TensorO1<double> phi1D_x(Np);  // Np Lagrange polynomials in 1D (x direction)
	TensorO1<double> phi1D_y(Np);  // Np Lagrange polynomials in 1D (y direction)

	TensorO1<double> phi_diff_1D_x(Np);  // Np derivatives of Lagrange polynomials in 1D
	TensorO1<double> phi_diff_1D_y(Np);  // Np derivatives of Lagrange polynomials in 1D

	this->LagrangePolynomial1D(x,xi,&phi1D_x); // Lagrange polynomials are computed in x dimension
	this->LagrangePolynomial1D(y,xi,&phi1D_y); // Lagrange polynomials are computed in y dimension

	this->LagrangePolynomialDerivative1D(x,xi,&phi_diff_1D_x); // Derivatives of Lagrange polynomials are computed in 1D in x
	this->LagrangePolynomialDerivative1D(y,xi,&phi_diff_1D_y); // Derivatives of Lagrange polynomials are computed in 1D in y

	// Partial derivative in x direction
	int counter = 0;
	for (Index i=0;i<Np;i++){
		for (Index j=0;j<Np;j++){
			phi_x->setValue(counter, phi_diff_1D_x.getValue(i)*phi1D_y.getValue(j)); //x: i loop, y: j loop
			counter++;
		};
	};

	// Partial derivative in y direction
	counter = 0;
	for (Index i=0;i<Np;i++){
		for (Index j=0;j<Np;j++){
			phi_y->setValue(counter, phi1D_x.getValue(i)*phi_diff_1D_y.getValue(j));
			counter++;
		};
	};
};


void FunctionalSpace::LagrangePolynomial3DTensor(double x, double y, double z, TensorO1<double> *xi, TensorO1<double> *phi){
	// (x,y,z): location at which interpolation is to be performed.      
	// (xi): An array of points where sampling of the data is done. Total number of elements of x = (N+1), where N = order of polynomial.
	// phi: array of the Lagrange polynomials phi[i] gives ith Lagrange polynomial

	
	int Np = xi->getSize();
	TensorO1<double> phi1Dx(Np); // Np 1D Lagrange polynomials in dim x
	TensorO1<double> phi1Dy(Np); // Np 1D Lagrange polynomials in dim y
	TensorO1<double> phi1Dz(Np); // Np 1D Lagrange polynomials in dim z
	this->LagrangePolynomial1D(x,xi,&phi1Dx); // Lagrange polynomials are computed for location x
	this->LagrangePolynomial1D(y,xi,&phi1Dy); // Lagrange polynomials are computed for location y
	this->LagrangePolynomial1D(z,xi,&phi1Dz); // Lagrange polynomials are computed for location z

	int counter = 0;
	for (Index i=0;i<Np;i++){
		for (Index j=0;j<Np;j++){
			for (Index k=0;k<Np;k++){
				double value = phi1Dx.getValue(i)*phi1Dy.getValue(j)*phi1Dz.getValue(k);
				phi->setValue(counter, value);
				counter++;
			};
		};
	};
};



void FunctionalSpace::LagrangePolynomialDerivative3DTensor(double x, double y, double z, TensorO1<double> *xi, TensorO1<double> *phi_x, TensorO1<double> *phi_y, TensorO1<double> *phi_z){
 	// The function fills derivatives of Lagrange Tensor polynomial in x, y and z direction
	// (x,y,z): location at which differentiation is to be performed.      
	// (xi): An array of points where sampling of the data is done. Total number of elements of x = N+1
	//    These points can be equi-spaced or roots of orthonormal polynomials. 
	// phi_x: array of the partial derivatives of Lagrange polynomials in x direction. phi_x[i] gives ith component of d(phi)/dx
	// phi_y: array of the partial derivatives of Lagrange polynomials in y direction. phi_y[i] gives ith component of d(phi)/dy
	// phi_z: array of the partial derivatives of Lagrange polynomials in z direction. phi_z[i] gives ith component of d(phi)/dz
	
	int Np = xi->getSize();


	// Lagrange polynomials in x, y and z directions
	TensorO1<double> phi1D_x(Np);  // Np Lagrange polynomials in x
	TensorO1<double> phi1D_y(Np);  // Np Lagrange polynomials in y
	TensorO1<double> phi1D_z(Np);  // Np Lagrange polynomials in z

	TensorO1<double> phi_DIFF_1D_x(Np); // Np derivatives of Lagrange polynomials in x
	TensorO1<double> phi_DIFF_1D_y(Np); // Np derivatives of Lagrange polynomials in y
	TensorO1<double> phi_DIFF_1D_z(Np); // Np derivatives of Lagrange polynomials in z

	this->LagrangePolynomial1D(x,xi,&phi1D_x); // Lagrange polynomials are computed in x dimension
	this->LagrangePolynomial1D(y,xi,&phi1D_y); // Lagrange polynomials are computed in y dimension
	this->LagrangePolynomial1D(z,xi,&phi1D_z); // Lagrange polynomials are computed in z dimension

	this->LagrangePolynomialDerivative1D(x,xi,&phi_DIFF_1D_x); // Derivatives of Lagrange polynomials are computed in x
	this->LagrangePolynomialDerivative1D(y,xi,&phi_DIFF_1D_y); // Derivatives of Lagrange polynomials are computed in y
	this->LagrangePolynomialDerivative1D(z,xi,&phi_DIFF_1D_z); // Derivatives of Lagrange polynomials are computed in z

	// Partial derivative in x direction
	int counter = 0;
	for (Index i=0;i<Np;i++){
		for (Index j=0;j<Np;j++){
			for (Index k=0;k<Np;k++){
				double value = phi_DIFF_1D_x.getValue(i) * phi1D_y.getValue(j) * phi1D_z.getValue(k);
				phi_x->setValue(counter, value);
				counter++;
			};
		};
	};

	// partial derivative in y direction
	counter = 0;
	for (Index i=0;i<Np;i++){
		for (Index j=0;j<Np;j++){
			for (Index k=0;k<Np;k++){
				double value = phi1D_x.getValue(i) * phi_DIFF_1D_y.getValue(j) * phi1D_z.getValue(k);
				phi_y->setValue(counter, value); 
				counter++;
			};
		};
	};

	// partial derivative in z direction
	counter = 0;
	for (Index i=0;i<Np;i++){
		for (Index j=0;j<Np;j++){
			for (Index k=0;k<Np;k++){
				double value =  phi1D_x.getValue(i) * phi1D_y.getValue(j) * phi_DIFF_1D_z.getValue(k);
				phi_z->setValue(counter, value); 
				counter++;
			};
		};
	};
};



void FunctionalSpace::generateLagMatrix1D(TensorO2<double> *LagMatrix){

	int nDOF = this->order +1;
	int nQuad = nDOF + this->intFlag; //number of quadrature points (=nDOF for inexact, Np+1 for exact)

	// Set the size of the Lagrange matrix (this re-allocates the memory)
	LagMatrix->setSize(nQuad, nDOF);

	// Generating DOF locations  (note: quadrature weights w_i not required. Only points x_i are used)
	TensorO1<double> x_i(nDOF);	// x_i are location of DOF points.
	TensorO1<double> w_i(nDOF);	// w_i are the quadrature weights. Will not be directly used.

	this->LGLRootsAndWeights1D(&x_i,&w_i); 

	// Generating quadrature points and weights
	TensorO1<double> x_q(nQuad);	// x_q are the location of quadrature points. 
	TensorO1<double> w_q(nQuad);	// w_q are the quadrature weights.

	this->LGLRootsAndWeights1D(&x_q, &w_q);  //for inexact integration, quadrature points are same as DOF locations.

	TensorO1<double> phi(nDOF);	//Temperory storage array for Lagrange polynomials at each quad point. 

	// Generating nDOF Lagrange polynomials for nQuad quadrature points (stored in nQuad X nDOF array)
	for (int k=0; k<nQuad; k++){
		this->LagrangePolynomial1D(x_q.getValue(k), &x_i, &phi);
		// now we have nDOF Lagrange polynomials based on x_i nodal points defined at the point x_q[k]
		for (int i=0; i<nDOF; i++){
			LagMatrix->setValue(k,i, phi.getValue(i));
		};
	};
};


// Lagrange polynomial matrix. first index runs over quadrature points, second over basis function numbers. 
// So for L, rows: quadrature points,  columns: Lagrange polynomials.
// For 2D implementation, quad points as well as weights are arranged in one straight row.
void FunctionalSpace::generateLagMatrix2DTensor(TensorO2<double>  *LagMatrix){
	// Basis functions computed using tensor products of 1D polynomials

	int Np = this->order +1;
	int K = Np + this->intFlag; //number of quadrature points (=Np for inexact, Np+1 for exact)

	int nDOF = Np * Np; //number of DOF points.
	int nQuad = K*K;    //number of quadrature points. For inexact integration, nQuad = nDOF

	// Set the size of the Lagrange matrix (this re-allocates the space for the matrix)
	LagMatrix->setSize(nQuad, nDOF);

	// Generating DOF locations on 1D element [-1,1] (note: quadrature weights w_i not required. Only points x_i are used)
	TensorO1<double> x_i(Np);
	TensorO1<double> w_i(Np);

	this->LGLRootsAndWeights1D(&x_i, &w_i); 


	// Generating quadrature points and weights on 2D standard element
	TensorO1<double> x_q(nQuad);	// x_q and y_q are location of quadrature points (x,y).
	TensorO1<double> y_q(nQuad);	
	TensorO1<double> w_q(nQuad);	// w_q are the quadrature weights.

	this->LGLRootsAndWeights2D(K, K, &x_q, &y_q, &w_q);  //for inexact integration, quadrature points are same as DOF locations.

	//Temperory storage array for Lagrange polynomials at each quad point.
	TensorO1<double> phi(nDOF);
	                        
	// Generating Np^2 Lagrange polynomials for K^2 quadrature points (stored in K^2 X Np^2 array)
	int quadPoint = 0;
	for (Index k1=0; k1<K; k1++){ // in x
		for (Index k2=0; k2<K; k2++){ //in y
			this->LagrangePolynomial2DTensor(x_q.getValue(quadPoint), y_q.getValue(quadPoint),  &x_i, &phi);
			// now we have Np*Np Lagrange polynomials based on x_i nodal points (in x and y) defined at the quadrature point
			for (Index DOFPoint=0; DOFPoint<Np*Np; DOFPoint++){
				LagMatrix->setValue(quadPoint, DOFPoint, phi.getValue(DOFPoint));
			};
			quadPoint ++;
		};
	};
};



//Lagrange polynomial matrix. first index runs over quadrature points, second over basis function numbers. 
//So for L, rows: quadrature points,  columns: Lagrange polynomials
//For 3D implementation, quad points as well as weights are arranged in one straight row.
void FunctionalSpace::generateLagMatrix3DTensor(TensorO2<double> *LagMatrix){

	int Np = this->order +1;
	int K = Np + this->intFlag; //number of quadrature points (=Np for inexact, Np+1 for exact)

	int nDOF = Np * Np * Np; //number of DOF points for 3D space.
	int nQuad = K*K*K;    //number of quadrature points for 3D space. For inexact integration, nQuad = nDOF.

	// Set the size of the Lagrange matrix (this re-allocates the memory)
	LagMatrix->setSize(nQuad, nDOF);

	// Generating DOF locations on 1D element [-1,1] (note: quadrature weights w_i not required. Only points x_i are used)
	TensorO1<double> x_i(Np);  //x_i locations of DOF on element [-1,1]. 
	TensorO1<double> w_i(Np);  //w_i will not be directly used (w_i = w_q for inexact integration). 

	this->LGLRootsAndWeights1D(&x_i, &w_i); 


	// Generating quadrature points and weights on 2D standard element
	TensorO1<double> x_q(nQuad);	// x coordinate of location of all quadrature points.
	TensorO1<double> y_q(nQuad);	// y coordinate of location of all quadrature points.
	TensorO1<double> z_q(nQuad);	// z coordinate of location of all quadrature points.
	TensorO1<double> w_q(nQuad);	// Quadrature weights arranged in a single array.
  	
	// For inexact integration, quadrature points are same as DOF locations. w_q will not be used here.
	this->LGLRootsAndWeights3D(K, K, K,  &x_q, &y_q, &z_q, &w_q);

 	// Temperory storage array for Lagrange polynomials at each quad point. 
	TensorO1<double> phi(nDOF);
	
	// Generating Np^3 Lagrange polynomials for K^3 quadrature points (stored in K^2 X Np^2 array)
	int quadPoint = 0;
	for (Index k1=0; k1<K; k1++){ // in x
		for (Index k2=0; k2<K; k2++){ //in y
			for (Index k3=0; k3<K; k3++){ //in z
				this->LagrangePolynomial3DTensor(x_q.getValue(quadPoint), y_q.getValue(quadPoint), z_q.getValue(quadPoint), &x_i, &phi);
				// now we have Np*Np*Np Lagrange polynomials based on x_i nodal points (in x, y and z) 

				for (int DOFPoint=0; DOFPoint<nDOF; DOFPoint++){
					LagMatrix->setValue(quadPoint, DOFPoint, phi.getValue(DOFPoint));
				};
				quadPoint ++;
			};
		};
	};
};




// Lagrange polynomial derivatives matrix. first index runs over quadrature points, second over derivatives of basis polynomial . 
// This stores values of derivatives of Lagrange polynomials at all quadrature points. The matrix has size nQuad X nDOF
// Rows: quadrature points,  columns: Lagrange polynomials derivatives
void FunctionalSpace::generateGradLagMatrix1D(TensorO2<double> *gradLagMatrix){

	int nDOF = this->order +1;
	int nQuad = nDOF + this->intFlag; //number of quadrature points (=nDOF for inexact, Np+1 for exact)

	// Set the size of the Lagrange matrix (this re-allocates the memory)
	gradLagMatrix->setSize(nQuad, nDOF);

	// Generating DOF locations  (note: quadrature weights w_i not required. Only points x_i are used)
	TensorO1<double> x_i(nDOF);	// x_i are location of DOF points.
	TensorO1<double> w_i(nDOF);	// w_i are the quadrature weights. Will not be directly used.

	this->LGLRootsAndWeights1D(&x_i,&w_i); 

	// Generating quadrature points and weights
	TensorO1<double> x_q(nQuad);	// x_q are the location of quadrature points. 
	TensorO1<double> w_q(nQuad);	// w_q are the quadrature weights.

	this->LGLRootsAndWeights1D(&x_q, &w_q);  //for inexact integration, quadrature points are same as DOF locations.

	TensorO1<double> phi_dx(nDOF);	//Temperory storage array for Lagrange polynomial derivatives at each quad point. 

	// Generating nDOF Lagrange polynomials for nQuad quadrature points (stored in nQuad X nDOF array)
	for (int k=0; k<nQuad; k++){
		this->LagrangePolynomialDerivative1D(x_q.getValue(k), &x_i, &phi_dx);
		// now we have nDOF Lagrange polynomials based on x_i nodal points defined at the point x_q[k]
		for (int i=0; i<nDOF; i++){
			gradLagMatrix->setValue(k,i, phi_dx.getValue(i));
		};
	};
};

void FunctionalSpace::generateGradLagMatrix2DTensor(TensorO2<double> *LagDrMatrix, TensorO2<double> *LagDsMatrix ){
	// Reference coordinate space: r,s space

	int Np = this->order +1;
	int K = Np + this->intFlag; //number of quadrature points (=Np for inexact, Np+1 for exact)

	int nDOF = Np * Np; //number of DOF points.
	int nQuad = K*K;    //number of quadrature points. For inexact integration, nQuad = nDOF

	// Set the size of the Lagrange matrix (this re-allocates the space for the matrix)
	LagDrMatrix->setSize(nQuad, nDOF);
	LagDsMatrix->setSize(nQuad, nDOF);

	// Generating DOF locations on 1D element [-1,1] (note: quadrature weights w_i not required. Only points x_i are used)
	TensorO1<double> x_i(Np);
	TensorO1<double> w_i(Np);

	this->LGLRootsAndWeights1D(&x_i, &w_i); 


	// Generating quadrature points and weights on 2D standard element
	TensorO1<double> x_q(nQuad);	// x_q and y_q are location of quadrature points (x,y).
	TensorO1<double> y_q(nQuad);	
	TensorO1<double> w_q(nQuad);	// w_q are the quadrature weights.

	this->LGLRootsAndWeights2D(K, K, &x_q, &y_q, &w_q);  //for inexact integration, quadrature points are same as DOF locations.

	//Temperory storage array for Lagrange polynomial derivatives at each quad point.
	TensorO1<double> phi_dr(nDOF);
	TensorO1<double> phi_ds(nDOF);
	                        
	// Generating Np^2 Lagrange polynomial derivatives for K^2 quadrature points (stored in K^2 X Np^2 arrays)
	int quadPoint = 0;
	for (Index k1=0; k1<K; k1++){ // in x
		for (Index k2=0; k2<K; k2++){ //in y
			this->LagrangePolynomialDerivative2DTensor(x_q.getValue(quadPoint), y_q.getValue(quadPoint),  &x_i, &phi_dr, &phi_ds);
			// now we have Np*Np Lagrange polynomials based on x_i nodal points (in x and y) defined at the quadrature point
			for (Index DOFPoint=0; DOFPoint<Np*Np; DOFPoint++){
				LagDrMatrix->setValue(quadPoint, DOFPoint, phi_dr.getValue(DOFPoint));
				LagDsMatrix->setValue(quadPoint, DOFPoint, phi_ds.getValue(DOFPoint));
			};
			quadPoint ++;
		};
	};
};

// Lagrange polynomial matrix. first index runs over quadrature points, second over basis function numbers. 
// So for L, rows: quadrature points,  columns: Lagrange polynomials
// For 3D implementation, quad points as well as weights are arranged in one straight row.
void FunctionalSpace::generateGradLagMatrix3DTensor(TensorO2<double> *LagDrMatrix, TensorO2<double> *LagDsMatrix, TensorO2<double> *LagDtMatrix){

	int Np = this->order +1;
	int K = Np + this->intFlag; //number of quadrature points (=Np for inexact, Np+1 for exact)

	int nDOF = Np * Np * Np; //number of DOF points for 3D space.
	int nQuad = K*K*K;    //number of quadrature points for 3D space. For inexact integration, nQuad = nDOF.

	// Set the size of the Lagrange matrix (this re-allocates the memory)
	LagDrMatrix->setSize(nQuad, nDOF);
	LagDsMatrix->setSize(nQuad, nDOF);
	LagDtMatrix->setSize(nQuad, nDOF);

	// Generating DOF locations on 1D element [-1,1] (note: quadrature weights w_i not required. Only points x_i are used)
	TensorO1<double> x_i(Np);  //x_i locations of DOF on element [-1,1]. 
	TensorO1<double> w_i(Np);  //w_i will not be directly used (w_i = w_q for inexact integration). 

	this->LGLRootsAndWeights1D(&x_i, &w_i); 


	// Generating quadrature points and weights on 2D standard element
	TensorO1<double> x_q(nQuad);	// x coordinate of location of all quadrature points.
	TensorO1<double> y_q(nQuad);	// y coordinate of location of all quadrature points.
	TensorO1<double> z_q(nQuad);	// z coordinate of location of all quadrature points.
	TensorO1<double> w_q(nQuad);	// Quadrature weights arranged in a single array.
  	
	// For inexact integration, quadrature points are same as DOF locations. w_q will not be used here.
	this->LGLRootsAndWeights3D(K, K, K,  &x_q, &y_q, &z_q, &w_q);

	//Temperory storage array for Lagrange polynomial derivatives at each quad point.
	TensorO1<double> phi_dr(nDOF);
	TensorO1<double> phi_ds(nDOF);
	TensorO1<double> phi_dt(nDOF);
	
	// Generating Np^3 Lagrange derivative polynomials for K^3 quadrature points (stored in K^2 X Np^2 array)
	int quadPoint = 0;
	for (Index k1=0; k1<K; k1++){ // in x
		for (Index k2=0; k2<K; k2++){ //in y
			for (Index k3=0; k3<K; k3++){ //in z
				this->LagrangePolynomialDerivative3DTensor(x_q.getValue(quadPoint), y_q.getValue(quadPoint), z_q.getValue(quadPoint), &x_i, &phi_dr, &phi_ds, &phi_dt);
				// now we have Np*Np*Np Lagrange polynomials based on x_i nodal points (in x, y and z) 

				for (int DOFPoint=0; DOFPoint<nDOF; DOFPoint++){
					LagDrMatrix->setValue(quadPoint, DOFPoint, phi_dr.getValue(DOFPoint));
					LagDsMatrix->setValue(quadPoint, DOFPoint, phi_ds.getValue(DOFPoint));
					LagDtMatrix->setValue(quadPoint, DOFPoint, phi_dt.getValue(DOFPoint));
				};
				quadPoint ++;
			};
		};
	};
};


//************************************************************************************************//
// Modal Space
//************************************************************************************************//


void FunctionalSpace::LegendrePolynomial1D(double r, int N, TensorO1<double> *phi){
	// The polynomials are defined on the reference element [-1,1]. 
	// Need to multiply by the Jacobian of transformation for the real element.
	// N: polynomial order
	// The polynomial value at location r and the first two derivatives are stored in phi
	// phi[0]: polynomial of order N at location r
	// phi[1]: first derivative;   phi[2]: second derivative;


	int Np = N + 1; //No of points = polynomial order + 1  for 1D



	// Polynomial at time level i (running index) 
	double L_0[3]; //0: polynomial, 1: first derivative, 2: second derivative
	// Polynomial at time level i - 1 (i.e. previous mode)
	double L_1[3]; //0: polynomial, 1: first derivative, 2: second derivative
	// Polynomial at time level i - 2 (i.e. previous to previous mode)
	double L_2[3]; //0: polynomial, 1: first derivative, 2: second derivative

	// initializing all array with zeros
	for (int i=0; i< 3; i++){
		L_0[i] = 0.0;
		L_1[i] = 0.0;
	};

	// setting the first polynomial to 1.0
	L_0[0] = 1.0;

	double a , b; //variables used in the calculation

	// recursive formula used.
	for (Index i=1; i< Np; i++){
		for (int k=0; k<3; k++){
			L_2[k] = L_1[k];
			L_1[k] = L_0[k];
		};

		a = (2.0* i - 1.0)/i;
		b = (i - 1.0)/ i;
		L_0[0] = a * r * L_1[0] - b * L_2[0];  // Recursive formula for Legendre polynomials. Refer class notes. 
		L_0[1] = a * (L_1[0] + r * L_1[1]) - b * L_2[1];  //differentiating polynomial formula once
		L_0[2] = a * (2.* L_1[1] + r * L_1[2]) - b * L_2[2];  //differentiating polynomial formula twice
	};

	for (Index k=0; k<3; k++){
		//phi[k] = L_0[k]; // Latest values copied
		
		// Normalization of the polynomials
		// Due to this, we get orthonormal polynomials
		double gamman = 2.0 / (2.0*N + 1.0);
		//double gamma = 1.0 / sqrt(gamman); //normalization function
		double gamma = 1.0; //no normalization function
		phi->setValue(k, L_0[k] * gamma); // Latest values copied and multiplied by the normalization function
	};
};


void FunctionalSpace::LegendrePolynomial3DTensor(double x, double y, double z, int Nx, int Ny, int Nz, TensorO1<double> *phi){
	// Jacobian of two dimensional partial differentials
	TensorO1<double> phi1Dx(3);
	TensorO1<double> phi1Dy(3);
	TensorO1<double> phi1Dz(3);
	
	this->LegendrePolynomial1D(x, Nx, &phi1Dx);
	this->LegendrePolynomial1D(y, Ny, &phi1Dy);
	this->LegendrePolynomial1D(z, Nz, &phi1Dz);
	
	int counter = 0;
	for (Index i=0;i<3;i++){
		for (Index j=0;j<3;j++){
			for (Index k=0;k<3;k++){
				phi->setValue(counter, phi1Dx.getValue(i)*phi1Dy.getValue(j)*phi1Dz.getValue(k));
				counter++;
			};
		};
	};
};





//************************************************************************************************//
// Quadrature roots and weights
//************************************************************************************************//

void FunctionalSpace::LGLRootsAndWeights1D(TensorO1<double> *x_lgl_roots, TensorO1<double> *w){
	// Computes roots of Legendre polynomials of order N-1 and adds two end points
	
	int Np = x_lgl_roots->getSize();
	assert (Np > 1 && "minimum two points required"); //minimum two points required
	int N = Np - 1;  // N: order of the polynomial, Np: number of points required
	int N_half = Np/2; //returns the floor value of the half domain
	
	double x, dx, x_next;
	TensorO1<double> phi(3);  // Legendre polynomials and first two derivatives

	for (Index i=1; i< N_half + 1; i++){
		// initialization using Chebyshev roots. Note: the Chebyshev roots in the notes are given for index starting from 0. 
		// here we start index from 1. thus the formula is adjusted for this change.
		x = cos((2.0*i - 1.0) * PI/(2.0* N + 2)); //confirm this formula once.  It could be 2.0*N + 1 as in the sample code.

		for (int k=1; k<21; k++){
			this->LegendrePolynomial1D(x, N, &phi);
			// LGL roots are -1 and 1 at the ends, and roots of N-1 Legendre polynomial in the middle
			dx = -(1 - x*x)*phi.getValue(1) / (-2.0 * x * phi.getValue(1) + (1-x*x)*phi.getValue(2)); //derivatives used.
			x = x + dx;
			if (sqrt(pow(x,2)) < NANO){
				break;
			};
		};
		
		x_lgl_roots->setValue(Np-i,  x);
		w->setValue(Np-i, 2.0 / (N*Np*phi.getValue(0)*phi.getValue(0) ));
	};

	// check for zero root
	if (N+1 != 2*N_half){
		x = 0.0;
		this->LegendrePolynomial1D(x, N, &phi);
		x_lgl_roots->setValue(N_half,  x);
		w->setValue(N_half,  2.0 / (N*Np*phi.getValue(0)*phi.getValue(0) ));
	};

	// Rest of the roots found using symmetry 
	for (Index i=1; i< N_half + 1; i++){
		x_lgl_roots->setValue(i-1,  -x_lgl_roots->getValue(Np-i)) ;
		w->setValue(i-1,            w->getValue(Np-i)) ;
	};
};



void FunctionalSpace::LGLRootsAndWeights2D(int Npx, int Npy, TensorO1<double> *x_lg_roots, TensorO1<double> *y_lg_roots, TensorO1<double> *w){


	TensorO1<double> xi_x(Npx);	// 1D roots for x dim
	TensorO1<double>  w_x(Npx);	//  Weights for x dim

	TensorO1<double> xi_y(Npy);	// 1D roots for y dim
	TensorO1<double>  w_y(Npy);	//  Weights for y dim

	this->LGLRootsAndWeights1D(&xi_x, &w_x);
	this->LGLRootsAndWeights1D(&xi_y, &w_y);
	
	int counter = 0;
	for (Index i=0;i<Npx;i++){
		for (Index j=0;j<Npy;j++){
			x_lg_roots->setValue(counter, xi_x.getValue(i));
			y_lg_roots->setValue(counter, xi_y.getValue(j));
			w->setValue(counter, w_x.getValue(i)*w_y.getValue(j));
			counter++;
		};
	};
};


void FunctionalSpace::LGLRootsAndWeights3D(int Npx, int Npy, int Npz, TensorO1<double> *x_lg_roots, TensorO1<double> *y_lg_roots, TensorO1<double> *z_lg_roots, TensorO1<double> *w){


	TensorO1<double> xi_x(Npx);	// 1D roots for x dim
	TensorO1<double>  w_x(Npx);	//  Weights for x dim

	TensorO1<double> xi_y(Npy);	// 1D roots for y dim
	TensorO1<double>  w_y(Npy);	//  Weights for y dim

	TensorO1<double> xi_z(Npz);	// 1D roots for z dim
	TensorO1<double>  w_z(Npz);	//  Weights for z dim

	this->LGLRootsAndWeights1D(&xi_x, &w_x);
	this->LGLRootsAndWeights1D(&xi_y, &w_y);
	this->LGLRootsAndWeights1D(&xi_z, &w_z);
	

	int counter = 0;
	for (Index i=0;i<Npx;i++){
		for (Index j=0;j<Npy;j++){
			for (Index k=0;k<Npz;k++){
				x_lg_roots->setValue(counter, xi_x.getValue(i));
				y_lg_roots->setValue(counter, xi_y.getValue(j));
				z_lg_roots->setValue(counter, xi_z.getValue(k));
				w->setValue(counter, w_x.getValue(i)*w_y.getValue(j)*w_z.getValue(k));

				counter++;
			};
		};
	};
};






//************************************************************************************************//
// Cell Matrices
//************************************************************************************************//

void FunctionalSpace::generateMassMatrix3DTensor(Cell* cell){
	
	int Np = this->order +1;
	int K = Np + this->intFlag; //number of quadrature points (=Np for inexact, Np+1 for exact)

	int nDOF = Np * Np * Np; //number of DOF points.
	assert(nDOF == cell->getNDOF());
	int nQuad = K*K*K;    //number of quadrature points. For inexact integration, nQuad = nDOF
	assert(nQuad == cell->getNQuad());

	// Computation of Mass matrix

	double quadSum;
	for (Index i=0; i<nDOF; i++){
		for (Index j=0; j<nDOF; j++){
			//i,j th location of M

			//quadrature performed over nQuad points
			quadSum = 0;
			for (Index k=0; k<nQuad; k++){
				double J = cell->getJacobian()->getValue(k);
				double w_q = cell->getReferenceCell()->getQuadWeights()->getValue(k);
				double phi_i = cell->getReferenceCell()->getLagMatrix()->getValue(k,i);
				double phi_j = cell->getReferenceCell()->getLagMatrix()->getValue(k,j);

				quadSum += phi_i * phi_j * w_q * J; 
			};
			cell->getMassMatrix()->setValue(i,j, quadSum);
		};
	};
};


void FunctionalSpace::generateVandermondeMatrix3DTensor(Cell *cell){
	int Np = this->order +1;
	int nDOF = Np * Np * Np; //number of DOF points.
	assert(nDOF == cell->getNDOF());

	TensorO1<double> x(nDOF);
	TensorO1<double> y(nDOF);
	TensorO1<double> z(nDOF);
	TensorO1<double> w(nDOF);

	// Legendre polynomials phi[0] and all the derivatives
	TensorO1<double> phi(27);

	this->LGLRootsAndWeights3D(Np,Np,Np, &x, &y, &z, &w);

	for (Index i=0; i<nDOF; i++){	
		int counter = 0;
		for (Index k=0; k<Np; k++){
			for (Index l=0; l<Np; l++){
				for (Index m=0; m<Np; m++){
					this->LegendrePolynomial3DTensor(x.getValue(i), y.getValue(i), z.getValue(i), k, l, m, &phi);
					// Legendre polynomial in the 0th position, rest are derivatives
					double phi0 = phi.getValue(0); 
					cell->getVandermondeMatrix()->setValue(i, counter, phi0);
					counter ++;
				};
			};
		};
	};
};


void FunctionalSpace::generateDiffMatrix3DTensor(Cell *cell){

	// Gradient of the Lagrange matrix computed on the reference coordinate system i.e. r,s,t system (given by LagDr(st)Matrix)
	// intFlag: 0:inexact integration, 1:exact integration
	// formFlag: 0:weak formulation, 1:strong formulation
	//	    Default value is taken as 0 i.e. weak form to be used for DG method
	// J: array of det(J) at quadrature points for cells

	// dRST_by_dXYZ is the inverse Jacobian matrix. It provides metric of transformation 
	// It is given by:
	// 	[ [dr/dx, dr/dy, dr/dz], 
	//	  [ds/dx, ds/dy, ds/dz], 
	//	  [dt/dx, dt/dy, dt/dz]]    matrix for nQuad points
	// dim: nQuad x 3 x 3


	//computing number of DOF location and no of quadrature points
	int Np = this->order +1;
	int K = Np + this->intFlag; //number of quadrature points (=Np for inexact, Np+1 for exact)

	int nDOF = Np * Np * Np; //number of DOF points.
	int nQuad = K*K*K;       //number of quadrature points. For inexact integration, nQuad = nDOF

	// Computation of Diff matrix
	double quadSum_x;
	double quadSum_y;
	double quadSum_z;
	
	// Change the value to 1 if strong formulation is required.
	int formFlag = 0; 
	assert(formFlag == 0 or formFlag == 1 && "Use 0:weak formulation, 1:strong formulation of the differentiation matrix D");

	for (int i=0; i<nDOF; i++){
		for (int j=0; j<nDOF; j++){
			//i,j th location of D

			// using The strong or the weak form:
			int firstIndex, secondIndex;
			firstIndex = formFlag * j + (1-formFlag)* i; //i.e. i if formFlag = 0, weak form; else j
			secondIndex = (1-formFlag)*j + formFlag * i; //i.e. j if formFlag = 1, strong form; else i

			// Assertion for weak formulation
			assert(firstIndex == i and secondIndex == j);

			// i.e. when formFlag = 0, (weak) then firstIndex = i, secondIndex = j;
			// i.e. the integration is \int (dphi_i phi_j) 
			//quadrature performed over K^2 points
			quadSum_x = 0;
			quadSum_y = 0;
			quadSum_z = 0;

			for (int k=0; k<nQuad; k++){
				//separating metric terms at kth quadrature point:
				double drBydx = cell->getInverseJacobian()->getValue(k, 0, 0);
				double drBydy = cell->getInverseJacobian()->getValue(k, 0, 1);
				double drBydz = cell->getInverseJacobian()->getValue(k, 0, 2);

				double dsBydx = cell->getInverseJacobian()->getValue(k, 1, 0);
				double dsBydy = cell->getInverseJacobian()->getValue(k, 1, 1);
				double dsBydz = cell->getInverseJacobian()->getValue(k, 1, 2);

				double dtBydx = cell->getInverseJacobian()->getValue(k, 2, 0);
				double dtBydy = cell->getInverseJacobian()->getValue(k, 2, 1);
				double dtBydz = cell->getInverseJacobian()->getValue(k, 2, 2);


				double phi_i_dr = cell->getReferenceCell()->getLagDrMatrix()->getValue(k,firstIndex);
				double phi_i_ds = cell->getReferenceCell()->getLagDsMatrix()->getValue(k,firstIndex);
				double phi_i_dt = cell->getReferenceCell()->getLagDtMatrix()->getValue(k,firstIndex);

				// partial derivative of Lagrange polynomial in x direction
				double phi_i_X = phi_i_dr * drBydx +  phi_i_ds * dsBydx + phi_i_dt * dtBydx;
				// partial derivative of Lagrange polynomial in y direction
				double phi_i_Y = phi_i_dr * drBydy +  phi_i_ds * dsBydy + phi_i_dt * dtBydy;
				// partial derivative of Lagrange polynomial in z direction
				double phi_i_Z = phi_i_dr * drBydz +  phi_i_ds * dsBydz + phi_i_dt * dtBydz;

				double phi_j = cell->getReferenceCell()->getLagMatrix()->getValue(k, secondIndex);

				double w_q = cell->getReferenceCell()->getQuadWeights()->getValue(k);

				double J = cell->getJacobian()->getValue(k);

				quadSum_x +=  phi_i_X * phi_j * w_q * J;
				quadSum_y +=  phi_i_Y * phi_j * w_q * J;
				quadSum_z +=  phi_i_Z * phi_j * w_q * J;
			};
			
			cell->getDxMatrix()->setValue(i, j, quadSum_x);
			cell->getDyMatrix()->setValue(i, j, quadSum_y);
			cell->getDzMatrix()->setValue(i, j, quadSum_z);
		};
	};
};


void FunctionalSpace::generateFluxMatrix3DTensor(Cell *cell){
	// Face* faces[6]: This indicates that faces includes 6 memebers each of the type pointer to class Face
	// 		   Other allowed forms are Face** faces and  Face* faces[] 
	// faceRelationshipArray: tells the relationship of the current cell with each of the 
	// 			  defining faces. refer Cell::faceRelationshipArray[] for details
	// w_q2D= Array of quadrature weights at all the quadrature points in the correct sequence for 2D surface.
	// LagMatrix: Matrix of Lagrange polynomials computed at quadrature points of reference cell. dim: nQuad x nDOF  (refer this->generateLagMatrix())
	// Basis functions computed using tensor products of 1D polynomials

	// intFlag: 0:inexact integration, 1:exact integration

	int Np = this->order +1;
	int K = Np + this->intFlag; //number of quadrature points (=Np for inexact, Np+1 for exact)

	int nDOF = Np * Np * Np; //number of DOF points.
	int nQuadFace = K*K;     //number of quadrature points. For inexact integration, nQuad = nDOF

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

	TensorO1<int>* FCArray; //stores mapOwnerQuadPoints or neighbourQuadPoints depending on faceRelationshipArray.
	for (int f=0; f<6; f++){
		Face *face = cell->getDefiningFace(f);
		RefFace* refFace = face->getReferenceFace();
		assert(refFace->getFaceType() == FaceType::Quad);
		// we are at face f
	
		// copy the mapQuadPoints array based on owner or neighbour cell
		if (cell->faceRelationshipArray.getValue(f) == 0){
			// 0 means owner
			FCArray = face->getOwnerQuadPointsArray();
		}
		else{
			//1 means neighbour
			FCArray = face->getNeighbourQuadPointsArray();
		};

		// Actual computation for each face
		for (Index i=0; i<cell->getNDOF(); i++){
			for (Index j=0; j<cell->getNDOF(); j++){
				// i,j location in Flux matrix

				// Initialize the flux matrix to zero
				cell->getFMatrix()->setValue(f, i, j, 0.0);

				for (Index k=0; k<nQuadFace; k++){
					//quadrature
					int globalQPlocation = FCArray->getValue(k); // refer Face::mapOwnerQuadPoints 
					// and Face::mapNeighbourQuadPoints for finding out more.
					double value =  cell->getReferenceCell()->getLagMatrix()->getValue(globalQPlocation, i) * 
						        cell->getReferenceCell()->getLagMatrix()->getValue(globalQPlocation, j) *  
							face->getReferenceFace()->getQuadWeights()->getValue(k) * 
							face->getJacobian()->getValue(k) ;

					cell->getFMatrix()->addValue(f,i,j, value);
				};
			};
		};
	};
};









//*****************************************************************************************************//
// Tests

namespace Test{

	double lagrangeInterpolation1D(double r, TensorO1<double> *x, TensorO1<double> *f){
		// Returns interpolated value of function f(x) at location r
		// The polynomial is constructed on the reference element [-1, 1]. Multiply by the Jacobian of transformation for the real element.
		// r: location at which interpolation is to be performed.    r E [-1,1]  
		// x: An array of points where sampling of the data is done. Total number of elements of x = N+1
		//    These points can be equi-spaced or roots of orthonormal polynomials. 
		// f: Value of the function f(x) at location given by array x

		int Np = x->getSize();

		double sum = 0;

		TensorO1<double> phi(Np);
		FunctionalSpace Fn(Np-1, IntFlag::exact);
		Fn.LagrangePolynomial1D(r,x,&phi); // Lagrange polynomials are computed for location r

		for (int i=0; i<Np; i++){
			sum = sum + f->getValue(i) * phi.getValue(i); //Lagrange interpolation formula
		};

		return sum;
	}

	double lagrangeDifferentiation1D(double r, TensorO1<double> *x, TensorO1<double> *f){
		// Returns interpolated value of derivative of function f(x) at location r. i.e. f'(x)
		// The differentiation polynomial is constructed on the reference element [-1, 1]. Multiply by the Jacobian of transformation for the real element.
		// r: location at which interpolation is to be performed.    r E [-1,1]  
		// x: An array of points where sampling of the data is done. Total number of elements of x = N+1
		//    These points can be equi-spaced or roots of orthonormal polynomials. 
		// f: Value of the function f(x) at location given by array x

		int Np = x->getSize();

		double sum = 0;

		TensorO1<double> phi(Np);
		FunctionalSpace Fn(Np-1, IntFlag::exact);
		Fn.LagrangePolynomialDerivative1D(r,x,&phi); // Derivatives of Lagrange polynomials are computed for location r

		for (int i=0; i<Np; i++){
			sum = sum + f->getValue(i) * phi.getValue(i); //Lagrange interpolation formula using derivatives of the polynomials
		};

		return sum;
	}


	double lagrangeInterpolation2D(double x, double y, TensorO1<double> *xi, TensorO1<double> *f){
		// Returns interpolated value of function f(x,y) at location (x,y)
		// The polynomial is constructed on the reference element [-1, 1]. Multiply by the Jacobian of transformation for the real element.
		// x,y: location at which interpolation is to be performed.    x,y E [-1,1]  
		// xi: An array of points where sampling of the data is done. Total number of elements of (N+1)(N+1)
		//    These points can be equi-spaced or roots of orthonormal polynomials. 
		//    Currently, only one dimensional array xi is required. The same array can be used in all the dimensions due to symmetry.
		//    So required size of xi is Np only.
		// f: Value of the function f(x,y) at location given by array xi

		int Np = xi->getSize();

		double sum = 0;

		TensorO1<double> phi(Np*Np); // Np Lagrange polynomials at location r
		FunctionalSpace Fn(Np-1, IntFlag::exact);
		Fn.LagrangePolynomial2DTensor(x,y,xi,&phi); // Lagrange polynomials are computed for location r

		for (Index i=0; i<Np*Np; i++){
			sum = sum + f->getValue(i) * phi.getValue(i); //Lagrange interpolation formula
		};
		return sum;
	};

	double lagrangeInterpolation3D(double x, double y, double z,  TensorO1<double> *xi, TensorO1<double> *f){
		// Returns interpolated value of function f(x,y,z) at location (x,y,z)
		// The polynomial is constructed on the reference element [-1, 1]. Multiply by the Jacobian of transformation for the real element.
		// x,y,z: location at which interpolation is to be performed.    x,y,z E [-1,1]  
		// N: order of the polynomial. 
		// xi: An array of points where sampling of the data is done. Total number of elements of (N+1)(N+1)(N+1)
		//    These points can be equi-spaced or roots of orthonormal polynomials. 
		// f: Value of the function f(xi,yi,zi) at location given by array xi

		int Np = xi->getSize();

		double sum = 0;

		TensorO1<double> phi(Np*Np*Np); // Np Lagrange polynomials at location r
		FunctionalSpace Fn(Np-1, IntFlag::exact);
		Fn.LagrangePolynomial3DTensor(x,y,z,xi,&phi); // Lagrange polynomials are computed for location r

		for (int i=0; i<Np*Np*Np; i++){
			sum = sum + f->getValue(i) * phi.getValue(i); //Lagrange interpolation formula
		};
		return sum;
	};




	// Generate Lagrange polynomials and export in a file
	void TestLagrangePolynomial(){

		cout << "\nRunning tests for 1D Lagrange polynomials...\n";
		cout << "First 15 polynomials are plotted.\n";

		for (Index order = 1 ; order < 16; order ++){

			// The following function generates the plotting script for automated plotting.
			string name = "Lagrange";
			generatePlottingScript(name, order);
			FunctionalSpace F(order, IntFlag::exact);
			int N = order;
			int Np = order + 1;
			double y[100]; //points which are plotted actually. Values of polynomials found out at these points.
			double Fn[100][30]; //up to order 30 polynomials found out. increase the size if more required.
			TensorO1<double> fn(30);     // local array of Lagrange polynomials
			double dy = 2./(99);
			for (int j=0; j< 100; j++){
				y[j] = -1.0 + j*dy;
			};

			TensorO1<double> x(Np); //x nodes where polynomials are 1 or zero
			TensorO1<double> w(Np); //corresponding weights for quadrature. Not required in this test.

			F.LGLRootsAndWeights1D(&x,&w);

			for (int j=0; j< 100; j++){
				F.LagrangePolynomial1D(y[j],x, fn); 
				for (int k=0; k<Np; k++){
					Fn[j][k] = fn.getValue(k);
				};
			};

			// Plotting Lagrange polynomials in a file
			ofstream file;
			file.open("tests/Lagrange_polynomials.dat");
			if (file.is_open()){
				for (int j=0; j< 100; j++){
					for (int n=0; n< Np; n++){
						file << Fn[j][n] << " ";
					};
					file << endl;
				};
			};
			
			// Use python script for plotting these polynomials. The script is present in the folder tests.
			string command;
			command = "python tests/plot_Lagrange_polynomials.py";
			const char *c = command.c_str();
			int i = system(c);

			// Remove the used files and scripts from tests folder (only keep the results)
			command = "rm tests/plot_Lagrange_polynomials.py";
			c = command.c_str();
			i = system(c);

			command = "rm tests/Lagrange_polynomials.dat";
			c = command.c_str();
			i = system(c);

			// Print progress bar
			double done = 1.0*order/15.0*100;
			printProgressBar(done, "Check out tests/ for results");
		};
		cout << "\nThe results will be generated inside tests/ folder. Kindly inspect the polynomials quality.\n";


		cout << "\nRunning tests for 2D Lagrange polynomials...\n";
		cout << "First 2 order-polynomials are plotted.\n";
		for (Index order = 1 ; order < 3; order ++){
			// The following function generates the plotting script for automated plotting.
			generatePlottingScript2D("Lagrange2", order);
			FunctionalSpace F(order, IntFlag::exact);
			int N = order;
			int Np = order + 1;
			double x[100][100]; //points which are plotted actually. Values of polynomials found out at these points.
			double y[100][100]; //points which are plotted actually. Values of polynomials found out at these points.
			double Fn[100][100][30*30]; //up to order 30 polynomials found out. increase the size if more required.
			TensorO1<double> fn(Np*Np); // Local array
			double dy = 2./(99);
			for (int i=0; i< 100; i++){
				for (int j=0; j< 100; j++){
					x[i][j] = -1.0 + i*dy;
					y[i][j] = -1.0 + j*dy;
				};
			};

			TensorO1<double> xi(Np); //x nodes where polynomials are 1 or zero
			TensorO1<double> wi(Np); //corresponding weights for quadrature. Not required in this test.

			F.LGLRootsAndWeights1D(&xi,&wi);

			ofstream file;
			file.open("tests/Lagrange2_polynomials.dat");
			if (file.is_open()){
				for (int i=0; i< 100; i++){
					for (int j=0; j< 100; j++){
						F.LagrangePolynomial2DTensor(x[i][j],y[i][j], &xi, &fn); 
						for (int k=0; k<Np*Np; k++){
							Fn[i][j][k] = fn.getValue(k);
							file << Fn[i][j][k] << " ";
						};
						file << endl;
					};
				};
			};

			// Use python script for plotting these polynomials. The script is present in the folder tests.
			string command;
			command = "python tests/plot_Lagrange2_polynomials.py";
			const char *c = command.c_str();
			int i = system(c);

			// Remove the used files and scripts from tests folder (only keep the results)
			command = "rm tests/plot_Lagrange2_polynomials.py";
			c = command.c_str();
			i = system(c);

			command = "rm tests/Lagrange2_polynomials.dat";
			c = command.c_str();
			i = system(c);

			// Print progress bar
			double done = 1.0*order/2.0*100;
			printProgressBar(done, "Check out tests/ for results");
		};
		cout << "\nThe results will be generated inside tests/ folder. Kindly inspect the polynomials quality.\n";

	};

	// This function solves Assignment 1 from Prof Shiva ME 766 course
	// An analytical function (cos(pi x/2)) is given between [-1,1]
	// Need to find interpolated value at sampling points using the number of
	// interpolating points or order
	void TestAssignment1(){
		cout << "\nRunning tests for accuracy of 1D Lagrange polynomials with LGL roots...\n";
		generatePlottingScriptAssignment1("Assignment1");

		int order = 64;
		int sampling_points = 101;

		// Header printing
			

		ofstream file1;
		file1.open("tests/Assignment1_interpolation.dat");
		ofstream file2;
		file2.open("tests/Assignment1_differentiation.dat");
		for (int N=1; N< order+1; N++){
			FunctionalSpace F(N, IntFlag::exact);
			int Np = N + 1;
			
			// Quad points and weights
			TensorO1<double> xi(Np); 
			TensorO1<double> wi(Np); 
			// Analytical Expression given for ui
			TensorO1<double> ui(Np);

			// This is where the interpolation function is used to find the value
			double x[sampling_points];
			double u[sampling_points];
			double u_analytical[sampling_points];
			double error[sampling_points];
			double rms_error = 0.0;
			double max_error = 0.0;

			
			for (int i=0; i<sampling_points; i++){
				x[i] = -1.0 + i*2.0/(sampling_points-1);
			}

			F.LGLRootsAndWeights1D(&xi,&wi);

			for (int i=0; i<Np; i++){
				ui.setValue(i, cos(PI*xi.getValue(i)/2.0));
			}
			

			// Main function 
			double sum=0;
			for (int i=0; i< sampling_points; i++){
				u[i] = lagrangeInterpolation1D(x[i], &xi, &ui); 
				u_analytical[i] = cos(PI*x[i]/2.0);
				error[i] = pow((u_analytical[i]-u[i]),2.0);
				rms_error+=error[i];
				error[i] = pow(error[i],0.5);
				if(error[i]>max_error)
					max_error = error[i];
				sum += error[i];
			};	
			file1 << N << " " << sum/sampling_points << " " << pow((rms_error/sampling_points),0.5) << " " << max_error  << endl;

			
			// Testing Differentiation of Lagrange polynomial
			// Main function 
			sum=0;
			for (int i=0; i< sampling_points; i++){
				u[i] = lagrangeDifferentiation1D(x[i], &xi, &ui); 
				u_analytical[i] = -3.1415/2.0*sin(3.1415*x[i]/2.0);
				error[i] = pow((u_analytical[i]-u[i]),2.0);
				rms_error+=error[i];
				error[i] = pow(error[i],0.5);
				if(error[i]>max_error)
					max_error = error[i];
				sum += error[i];
			};	

			file2 << N << " " << sum/sampling_points << " " << pow((rms_error/sampling_points),0.5) << " " << max_error  << endl;

			printProgressBar(N/64.0*100, "Check out tests/ for results");
		};
		cout << endl;
		file1.close();
		file2.close();

		string command;
		command = "python tests/plot_error_Assignment1.py";
		const char *c = command.c_str();
		int i = system(c);

		// Remove the used files and scripts from tests folder (only keep the results)
		command = "rm tests/plot_error_Assignment1.py";
		c = command.c_str();
		i = system(c);

		assert(fileExists("tests/Assignment1_interpolation.dat"));
		command = "rm tests/Assignment1_interpolation.dat";
		c = command.c_str();
		i = system(c);

		assert(fileExists("tests/Assignment1_differentiation.dat"));
		command = "rm tests/Assignment1_differentiation.dat";
		c = command.c_str();
		i = system(c);

		cout << "\nThe results will be generated inside tests/ folder. Kindly inspect the error plots 1. interpolation_error.eps 2. differentiation_error.eps\n";


		//**************************************************************************************************//


		cout << "\nRunning tests for accuracy of 2D Lagrange polynomials with LGL roots...\n";
		generatePlottingScriptAssignment2("Assignment2");

		order = 40;
		sampling_points = 51;

		ofstream file3;
		file3.open("tests/Assignment2_interpolation.dat");
		for (int N=1; N< order+1; N++){
			FunctionalSpace F(N, IntFlag::exact);
			int Np = N + 1;
			
			// Quad points and weights
			TensorO1<double> xi(Np); 
			TensorO1<double> wi(Np); 
			// Analytical Expression given for ui
			TensorO1<double> ui(Np*Np);

			// This is where the interpolation function is used to find the value
			double x[sampling_points][sampling_points];
			double y[sampling_points][sampling_points];
			double u[sampling_points][sampling_points];
			double u_analytical[sampling_points][sampling_points];
			double error[sampling_points][sampling_points];      
			double rms_error = 0.0;
			double max_error = 0.0;

			
			for (int i=0; i<sampling_points; i++){
				for (int j=0; j<sampling_points; j++){
					x[i][j] = -1.0 + i*2.0/(sampling_points-1);
					y[i][j] = -1.0 + j*2.0/(sampling_points-1);
				}
			}
			

			F.LGLRootsAndWeights1D(&xi,&wi);


			int counter = 0;
			for (Index i=0; i<Np; i++){
				for (Index j=0; j<Np; j++){
					ui.setValue(counter, (1.0/pow(2*PI,0.5))*exp(-xi.getValue(i)*xi.getValue(i) - xi.getValue(j)*xi.getValue(j)));
					counter++;
				}
			}
			
			// Main function 
			double sum=0;
			for (Index i=0; i< sampling_points; i++){
				for (Index j=0; j< sampling_points; j++){
					u[i][j] = lagrangeInterpolation2D(x[i][j],y[i][j], &xi, &ui); 
					u_analytical[i][j] = (1.0/pow(2*PI,0.5))*exp(-x[i][j]*x[i][j] - y[i][j]*y[i][j]);
					error[i][j] = pow((u_analytical[i][j]-u[i][j]),2.0);
					rms_error+=error[i][j];
					error[i][j] = pow(error[i][j],0.5);
					if(error[i][j]>max_error)
						max_error = error[i][j];
					sum += error[i][j];
				};	
			};	


			file3 << N << " " << sum/sampling_points << " " << pow((rms_error/sampling_points),0.5) << " " << max_error  << endl;
			printProgressBar(N/40.0*100, "Check out tests/ for results");
		};
		cout << endl;
		file3.close();

		command = "python tests/plot_error_Assignment2.py";
		c = command.c_str();
		i = system(c);

		// Remove the used files and scripts from tests folder (only keep the results)
		command = "rm tests/plot_error_Assignment2.py";
		c = command.c_str();
		i = system(c);

		assert(fileExists("tests/Assignment2_interpolation.dat"));
		command = "rm tests/Assignment2_interpolation.dat";
		c = command.c_str();
		i = system(c);

		cout << "\nThe results will be generated inside tests/ folder. Kindly inspect the error plot: interpolation_error_Assignment2.eps\n";


		//**************************************************************************************************//


		cout << "\nRunning tests for accuracy of 3D Lagrange polynomials with LGL roots...\n";
		generatePlottingScriptAssignment2("Assignment3");

		order = 40;
		sampling_points = 10;

		ofstream file4;
		file4.open("tests/Assignment3_interpolation.dat");
		for (int N=1; N< order+1; N++){
			FunctionalSpace F(N, IntFlag::exact);
			int Np = N + 1;
			
			// Quad points and weights
			TensorO1<double> xi(Np); 
			TensorO1<double> wi(Np); 
			// Analytical Expression given for ui
			TensorO1<double> ui(Np*Np*Np);

			// This is where the interpolation function is used to find the value
			double x[sampling_points][sampling_points][sampling_points];
			double y[sampling_points][sampling_points][sampling_points];
			double z[sampling_points][sampling_points][sampling_points];
			double u[sampling_points][sampling_points][sampling_points];
			double u_analytical[sampling_points][sampling_points][sampling_points];
			double error[sampling_points][sampling_points][sampling_points];      
			double rms_error = 0.0;
			double max_error = 0.0;

			
			for (int i=0; i<sampling_points; i++){
				for (int j=0; j<sampling_points; j++){
					for (int k=0; k<sampling_points; k++){
						x[i][j][k] = -1.0 + i*2.0/(sampling_points-1);
						y[i][j][k] = -1.0 + j*2.0/(sampling_points-1);
						z[i][j][k] = -1.0 + k*2.0/(sampling_points-1);
					};
				}
			}
			

			F.LGLRootsAndWeights1D(&xi,&wi);

			int counter = 0;
			for (Index i=0; i<Np; i++){
				for (Index j=0; j<Np; j++){
					for (Index k=0; k<Np; k++){
						ui.setValue(counter, (1.0/pow(2*PI,0.5))*exp(-xi.getValue(i)*xi.getValue(i) - xi.getValue(j)*xi.getValue(j) - xi.getValue(k)));
						counter++;
					}
				}
			}

			
			// Main function 
			double sum=0;
			for (int i=0; i< sampling_points; i++){
				for (int j=0; j< sampling_points; j++){
					for (int k=0; k< sampling_points; k++){
						u[i][j][k] = lagrangeInterpolation3D(x[i][j][k],y[i][j][k],z[i][j][k],&xi, &ui); 
						u_analytical[i][j][k] = (1.0/pow(2*PI,0.5))*exp(-x[i][j][k]*x[i][j][k] - y[i][j][k]*y[i][j][k] - z[i][j][k]);
						error[i][j][k] = pow((u_analytical[i][j][k]-u[i][j][k]),2.0);
						rms_error+=error[i][j][k];
						error[i][j][k] = pow(error[i][j][k],0.5);
						if(error[i][j][k]>max_error)
							max_error = error[i][j][k];
						sum += error[i][j][k];
					};
				};	
			};	


			file4 << N << " " << sum/sampling_points << " " << pow((rms_error/sampling_points),0.5) << " " << max_error  << endl;
			printProgressBar(N/40.0*100, "Check out tests/ for results");
		};
		cout << endl;
		file4.close();

		command = "python tests/plot_error_Assignment3.py";
		c = command.c_str();
		i = system(c);

		// Remove the used files and scripts from tests folder (only keep the results)
		command = "rm tests/plot_error_Assignment3.py";
		c = command.c_str();
		i = system(c);

		assert(fileExists("tests/Assignment3_interpolation.dat"));
		command = "rm tests/Assignment3_interpolation.dat";
		c = command.c_str();
		i = system(c);

		cout << "\nThe results will be generated inside tests/ folder. Kindly inspect the error plot: interpolation_error_Assignment3.eps\n";


		//**************************************************************************************************//
	};

	void TestLGLRootsAndWeights(){

		cout << "\nRunning tests for 1D LGL Roots and weights...\n";

		ofstream file1;
		file1.open("tests/1D_LGL_Roots.dat");


		for (int N=1; N<20; N++){
			FunctionalSpace F(N, IntFlag::exact);
			double total_w = 0.0;
			int Np = N+1;
			file1 << "\n______________________________________________" << endl;
			file1 << "N: " << N << " Np: " << Np << endl;
			TensorO1<double> x(Np);
			TensorO1<double> w(Np);
			//double x[Np], w[Np];
			
			F.LGLRootsAndWeights1D(&x,&w);

			file1 << "roots \\in [-1,1]: " << endl;
			for (int n=0; n< Np; n++){
				file1 << x.getValue(n) << " ";
			};
			file1 << endl;
			file1 << "weights: " << endl;
			for (int n=0; n< Np; n++){
				total_w += w.getValue(n);
				file1 << w.getValue(n) << " ";
			};
			file1 << endl;
			file1<<"Total Weight : "<<total_w<<endl;
		};

		file1.close();

		cout << "1D LGL Roots and weights are printed in tests/1D_LGL_Roots.dat for inspection.\n";

		//************************************************************************************************//

		cout << "\nRunning tests for 2D LGL Roots and weights...\n";

		ofstream file2;
		file2.open("tests/2D_LGL_Roots.dat");

		for (int N=1; N<20; N++){
			FunctionalSpace F(N, IntFlag::exact);
			int Np = N+1;
			file2 << "\n______________________________________________" << endl;
			file2 << "N: " << N << " Np: " << Np << endl;
			TensorO1<double> x(Np*Np);
			TensorO1<double> y(Np*Np);
			TensorO1<double> w(Np*Np);
			
			F.LGLRootsAndWeights2D(Np,Np,&x,&y,&w);
			double total_w = 0.0;

			int counter = 0;
			for (int n=0; n< Np; n++){	
				for (int m=0; m< Np; m++){
					file2 << "("<<x.getValue(counter) << " "<<y.getValue(counter) << "), "<<w.getValue(counter)<<endl;
					total_w += w.getValue(counter);
					counter++;
				};
			};
			file2<<"Total Weight : "<<total_w<<endl;
		};
		file2.close();
		cout << "2D LGL Roots and weights are printed in tests/2D_LGL_Roots.dat for inspection.\n";

		//************************************************************************************************//

		cout << "\nRunning tests for 3D LGL Roots and weights...\n";

		ofstream file3;
		file3.open("tests/3D_LGL_Roots.dat");

		for (int N=1; N<20; N++){
			FunctionalSpace F(N, IntFlag::exact);
			int Np = N+1;
			file3 << "\n______________________________________________" << endl;
			file3 << "N: " << N << " Np: " << Np << endl;
			TensorO1<double> x(Np*Np*Np);
			TensorO1<double> y(Np*Np*Np);
			TensorO1<double> z(Np*Np*Np);
			TensorO1<double> w(Np*Np*Np);
			
			F.LGLRootsAndWeights3D(Np,Np,Np,&x,&y,&z,&w);
			double total_w = 0.0;

			int counter = 0;
			for (int n=0; n< Np; n++){	
				for (int m=0; m< Np; m++){
					for (int l=0; l< Np; l++){
						file3 << "("<<x.getValue(counter) << " "<<y.getValue(counter) << " "<<z.getValue(counter)  << "), "<<w.getValue(counter)<<endl;
						total_w += w.getValue(counter);
						counter++;
					};
				};
			};
			file3<<"Total Weight : "<<total_w<<endl;
		};
		file3.close();
		cout << "3D LGL Roots and weights are printed in tests/3D_LGL_Roots.dat for inspection.\n";
	};

	void TestLagrangeMatrix(){
		cout << "\nRunning tests for 1D Lagrange matrix...\n";
		cout << "1) Testing 'inexact' integration...\n";

		ofstream fileInexact1D;
		fileInexact1D.open("tests/LagrangeMatrix1DInexact.dat");
		for (Index order=1; order <= 10; order ++){
			FunctionalSpace F(order, IntFlag::inexact);
			
			int nDOF = order + 1; 
			int nQuad = order + 1 + F.getIntFlag();
			TensorO2<double> LagMatrix(nQuad, nDOF);

			F.generateLagMatrix1D(&LagMatrix);
			

			fileInexact1D << "Order: " << F.getOrder() << endl;
			for (Index i=0; i<nQuad; i++){
				for (Index j=0; j<nDOF; j++){
					fileInexact1D << fixed << setprecision(2) << LagMatrix.getValue(i,j) << "\t";
				};
				fileInexact1D << endl;
			};
			fileInexact1D << endl;
		};

		fileInexact1D.close();
		cout << "Matrices up to order 10 are printed inside tests/LagrangeMatrix1DInexact.dat. Kindly inspect the matrices.\n";
		cout << "2) Testing 'exact' integration...\n";

		ofstream fileExact1D;
		fileExact1D.open("tests/LagrangeMatrix1DExact.dat");
		for (Index order=1; order <= 10; order ++){
			FunctionalSpace F(order, IntFlag::exact);
			
			int nDOF = order + 1; 
			int nQuad = order + 1 + F.getIntFlag();
			TensorO2<double> LagMatrix(nQuad, nDOF);

			F.generateLagMatrix1D(&LagMatrix);
			

			fileExact1D << "Order: " << F.getOrder() << endl;
			for (Index i=0; i<nQuad; i++){
				for (Index j=0; j<nDOF; j++){
					fileExact1D <<  fixed << setprecision(2) << LagMatrix.getValue(i,j) << "\t";
				};
				fileExact1D << endl;
			};
			fileExact1D << endl;
		};

		fileExact1D.close();

		cout << "Matrices up to order 10 are printed inside tests/LagrangeMatrix1DExact.dat. Kindly inspect the matrices.\n";

		cout << "Test::TestLagrangeMatrix() passed.\n";
	};
};
