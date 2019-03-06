#ifndef REFCELL_H
#define REFCELL_H

#include "sysInclude.h"
#include "Point.h"
#include "FunctionalSpaces.h"

class RefCell{
	public:
		
		/// Constructor
		RefCell(); 

		/// Destructor
		~RefCell();

		/// Initialization
		void init(int order, IntFlag::intflag intFlag, CellType::cellType celltype) ;

		/// Set the integration type (exact or inexact)
		void setIntegrationType(IntFlag::intflag intFlag);

		/// Get the celltype 
		CellType::cellType getCellType();

		// Get and set methods for variables 
		int getNoOfPoints()const;
		void setNoOfPoints(int pts);

		int getNoOfFaces()const;
		void setNoOfFaces(int fcs);

		int getNoOfDOFPoints()const;
		void setNoOfDOFPoints(int dofs);

		/// Generate weights for Gaussian quadrature
		void generateQuadratureWeights();
		/// Generate matrix of Lagrange interpolation polynomials. Rows indicate varying points (DOF or Quad), columns indicate varying Lagrange polynomials
		void generateLagMatrix();
		/// Generate Gradient of Lagrange matrix
		void generateGradLagMatrix();
		//add generator functions for other matrices

		/// Get Quadrature weights
		TensorO1<double>* getQuadWeights();

		/// Get Lagrange matrix
		TensorO2<double>* getLagMatrix();

		/// Get Lagrange matrix
		TensorO2<double>* getLagDrMatrix();

		/// Get Lagrange matrix
		TensorO2<double>* getLagDsMatrix();

		/// Get Lagrange matrix
		TensorO2<double>* getLagDtMatrix();


		void printLagMatrix();
		void printGradLagMatrix();
		void print ();


	private:

		int order;	// Order of polynomial reconstruction (not order of accuracy of the resulting scheme)
		int noOfPoints; // No. of defining points
		int nDOF;       // No. of degrees of freedom. 
		int nQuad;      // No. of Quadrature points.
		int noOfFaces;  // No. of defining faces

		CellType::cellType celltype;

		IntFlag::intflag intFlag; //0:inexact, 1: exact;

		TensorO1<double> J;   //array of det(Jacobian)
		TensorO1<double> w_q; // array of quadrature weights at all the quadrature points.

		TensorO2<double> LagMatrix;    //Matrix of Lagrange polynomials. Size: nQuad X nDOF  (quadrature points i, polynomials j)
		TensorO2<double> LagDrMatrix ; //Derivative matrix of Lagrange polynomials in r direction (analogous to x in ref. plane)
		TensorO2<double> LagDsMatrix ; //Derivative matrix of Lagrange polynomials in s direction (analogous to y in ref. plane)
		TensorO2<double> LagDtMatrix ; //Derivative matrix of Lagrange polynomials in t direction (analogous to z in ref. plane)
};



namespace Test{
	void TestReferenceCells();
};

#endif
