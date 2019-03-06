#ifndef FUNCTIONALSPACES_H
#define FUNCTIONALSPACES_H

#include "sysInclude.h"
#include "mathFunctions.h"
#include "RefCell.h"
#include "RefFace.h"
#include "Face.h"
#include "Cell.h"

class Cell;


class FunctionalSpace{
	public:
		/// Constructor
		FunctionalSpace(int ord, IntFlag::intflag intFlag);

		/// Destructor
		~FunctionalSpace();

		/// Get order
		unsigned int getOrder()const;

		/// Set order
		void setOrder(unsigned int order);

		/// Set integration flag
		void setIntFlag(IntFlag::intflag intFlag);

		/// Get integration flag
		IntFlag::intflag getIntFlag();

		

		//********************************************************************************************************//
		// Functions related to nodal space
		//********************************************************************************************************//

		// The function fills the array of the N+1 Lagrange polynomials phi.
		// N is order of the polynomials. x is the array of locations of the sampling points.
		// r is the location where polynomials are required. 
		/// 1D Lagrange interpolation polynomial.
		void LagrangePolynomial1D(double r, TensorO1<double> *x, TensorO1<double> *phi); 	// Pass by pointer
		void LagrangePolynomial1D(double r, TensorO1<double> &x, TensorO1<double> &phi);  	// Pass by reference

 		// The function fills the array of derivatives of the N+1 Lagrange polynomials phi. i.e. phi'
		// r is the location where polynomials are required. 
		/// Derivative of 1D Lagrange polynomials
		void LagrangePolynomialDerivative1D(double r, TensorO1<double> &x, TensorO1<double> &phi); 	// Pass by reference
		void LagrangePolynomialDerivative1D(double r, TensorO1<double> *x, TensorO1<double> *phi); 	// Pass by reference

		// This function fills the array of Lagrange polynomials in a Tensor-product space.
		// This polynomials are required for quad-integration on the faces. (Or for a 2D quad based DG scheme)
		/// 2D Lagrange interpolation polynomials.
		void LagrangePolynomial2DTensor(double x, double y, TensorO1<double> *xi, TensorO1<double> *phi); 	// Pass by pointer

 		// This function fills two arrays of derivatives of the N+1 Lagrange polynomials phi. i.e. phi'_x and phi'_y
		// xi is the array of locations of the sampling points.
		// (x, y) is the location where polynomials are required. 
		/// Derivatives of 2D Lagrange polynomials
		void LagrangePolynomialDerivative2DTensor(double x, double y, TensorO1<double> *xi, TensorO1<double> *phi_x, TensorO1<double> *phi_y); 

		/// 3D Lagrange interpolation polynomials.
		void LagrangePolynomial3DTensor(double x, double y, double z, TensorO1<double> *xi, TensorO1<double> *phi);
		
		void LagrangePolynomialDerivative3DTensor(double x, double y, double z, TensorO1<double> *xi, TensorO1<double> *phi_x, TensorO1<double> *phi_y, TensorO1<double> *phi_z);


		/// Lagrange polynomial matrix 
		/// Stores pre-computed Lagrange polynomial values (all polynomials for the given order) at DOF (or Quadrature) locations.
		/// Dimension: nQuad x nDOF (i.e. rows indicate the quadrature points, columns indicate polynomials)
		/// LagMatrix->getValue(i,j) will return jth Lagrange polynomial interpolated at ith Quadrature point.
		void generateLagMatrix1D(TensorO2<double> *LagMatrix);

		/// Lagrange matrix for 2D Tensor product space
		void generateLagMatrix2DTensor(TensorO2<double>  *LagMatrix);

		/// Lagrange matrix for 3D Tensor product space
		void generateLagMatrix3DTensor(TensorO2<double>  *LagMatrix);

		/// Lagrange derivative matrix for 1D space
		void generateGradLagMatrix1D(TensorO2<double> *gradLagMatrix);

		/// Lagrange derivative matrix for 2D Tensor space
		void generateGradLagMatrix2DTensor(TensorO2<double> *LagDrMatrix, TensorO2<double> *LagDsMatrix );

		/// Lagrange derivative matrix for 3D Tensor space
		void generateGradLagMatrix3DTensor(TensorO2<double> *LagDrMatrix, TensorO2<double> *LagDsMatrix, TensorO2<double> *LagDtMatrix);

		//********************************************************************************************************//
		// Functions related to modal space
		//********************************************************************************************************//

		/// Returns Legendre polynomials and the first two derivatives at location $r \in [-1,1]$ for order N
		void LegendrePolynomial1D(double r, int N, TensorO1<double> *phi);

		/// 3D Legendre polynomial and all the derivatives. phi.getValue(0) returns the Legendre polynomial in 3D.
		void LegendrePolynomial3DTensor(double x, double y, double z, int Nx, int Ny, int Nz, TensorO1<double> *phi);

		//********************************************************************************************************//
		// Functions related to quadrature roots and weights
		//********************************************************************************************************//

		// Collocation points (roots of orthogonal modes) and quadrature weights
		// x : points array (of size Np)
		// w : weights
		/// Finds Legendre-Gauss_Lobatto roots and weights for 1D interpolation.
		void LGLRootsAndWeights1D(TensorO1<double> *x_lgl_roots, TensorO1<double> *w);

		//  This function implicitly uses 1D routine and takes a Tensor product for finding roots and weights for 2D.
		/// Finds LGL roots and weights for 2D interpolation.
		void LGLRootsAndWeights2D(int Npx, int Npy, TensorO1<double> *x_lg_roots, TensorO1<double> *y_lg_roots, TensorO1<double> *w);

		/// Finds LGL roots and weights for 3D interpolation.
		void LGLRootsAndWeights3D(int Npx, int Npy, int Npz, TensorO1<double> *x_lg_roots, TensorO1<double> *y_lg_roots, TensorO1<double> *z_lg_roots, TensorO1<double> *w);



		//********************************************************************************************************//
		// Functions related to Cell Matrices
		//********************************************************************************************************//

		/// Generate Mass matrix for Hex cell 
		void generateMassMatrix3DTensor(Cell* cell);
		/// Generate Vandermonde matrix for Hex cell 
		void generateVandermondeMatrix3DTensor(Cell *cell);
		/// Generate Differentiation matrix for Hex cell 
		void generateDiffMatrix3DTensor(Cell *cell);
		/// Generate Flux matrix for Hex cell 
		void generateFluxMatrix3DTensor(Cell *cell);





	private:

		/// Order is the order of the interpolation not the numerical order of accuracy of the scheme.
		/// i.e. if order == 0: piecewise constant,  order == 1: linear etc.
		unsigned int order;

		/// Integration flag, exact or inexact integration. Refer intFlag.h for details.
		IntFlag::intflag intFlag;
};



namespace Test{

		double lagrangeInterpolation1D(double r, TensorO1<double> *x, TensorO1<double> *f); // interpolates function f(x) sampled at x points (i.e. array f at points x) at location r using Nth order Lagrange polynomial.
		double lagrangeDifferentiation1D(double r, TensorO1<double> *x, TensorO1<double> *f); // interpolates derivatives of function f(x) sampled at x points (i.e. array f at points x) at location r using Nth order Lagrange polynomial.

		double lagrangeInterpolation2D(double x, double y, TensorO1<double> *xi, TensorO1<double> *f);
		double lagrangeInterpolation3D(double x, double y, double z, TensorO1<double> *xi, TensorO1<double> *f);

		void TestLagrangePolynomial();
		void TestAssignment1();
		void TestLGLRootsAndWeights();

		void TestLagrangeMatrix();
};
#endif
