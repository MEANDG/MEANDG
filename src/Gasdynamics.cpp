#include "Gasdynamics.h"


void getConservedVariableVector(TensorO1<double> *primitiveVariableVector, TensorO1<double> *conservedVariableVector){
	//primitive variable vector needs to be in the following order:
	//[rho, u, v, w, p];
	//where, rho: density;    {u,v,w}: velocity components in cartesian x,y,z coord;      p: pressure

	double rho, u, v, w, p;
	//extract variables
	rho = primitiveVariableVector->getValue(0);
	if (abs(rho) < EPS){
		rho = abs(rho) + EPS;
	};
	u = primitiveVariableVector->getValue(1);
	v = primitiveVariableVector->getValue(2);
	w = primitiveVariableVector->getValue(3);
	p = primitiveVariableVector->getValue(4);


	// construct energy:

	double e;  //note that e is actually total energy (rho*E, where E is specific total energy) 
		   //so no need to multiply by rho for conserved variable vector


	e = p/(GAMMA -1 )  + 1.0/2.0 * rho * (u*u + v*v + w*w);

	//assemble conserved variable vector:
	conservedVariableVector->setValue(0, rho);
	conservedVariableVector->setValue(1, rho * u);
	conservedVariableVector->setValue(2, rho * v);
	conservedVariableVector->setValue(3, rho * w);
	conservedVariableVector->setValue(4, e);
};


void getPrimitiveVariableVector(TensorO1<double> *conservedVariableVector, TensorO1<double> *primitiveVariableVector){
	//conserved variable vector needs to be in the following order:
	//[rho, rhou, rhov, rhow, e]

	double rho, u, v, w, p, e;

	// extract variables
	rho = conservedVariableVector->getValue(0); //so that when we divide by rho, and if rho=0, then there shouldn't be a problem.
	if (abs(rho) < EPS){
		rho = abs(rho) + EPS;
	};

	u   = conservedVariableVector->getValue(1) / (rho); 
	v   = conservedVariableVector->getValue(2) / (rho); 
	w   = conservedVariableVector->getValue(3) / (rho); 
	e   = conservedVariableVector->getValue(4);
	p = (GAMMA - 1) * (e  -   (1.0/2.0 * rho *(u*u + v*v + w*w) ) );


	//assemble primitive variable vector:
	primitiveVariableVector->setValue(0, rho);
	primitiveVariableVector->setValue(1, u);
	primitiveVariableVector->setValue(2, v);
	primitiveVariableVector->setValue(3, w);
	primitiveVariableVector->setValue(4, p);
};



void getEulerFluxFromConservedVariables(TensorO1<double> *conservedVariableVector, TensorO1<double> *fluxVectorX, TensorO1<double> *fluxVectorY, TensorO1<double> *fluxVectorZ ){
	double rho, u, v, w, p, e;
	double GAMMA = 1.4;

	// extract variables
	rho = conservedVariableVector->getValue(0); //so that when we divide by rho, and if rho=0, then there shouldn't be a problem.
	u   = conservedVariableVector->getValue(1) / (rho+NANO); 
	v   = conservedVariableVector->getValue(2) / (rho+NANO); 
	w   = conservedVariableVector->getValue(3) / (rho+NANO); 
	e   = conservedVariableVector->getValue(4);
	p = (GAMMA - 1) * (e  -   (1.0/2.0 * rho *(u*u + v*v + w*w) ) );


	// assemble flux vector in x direction:
	fluxVectorX->setValue(0, rho * u);
	fluxVectorX->setValue(1, rho * u*u  + p);
	fluxVectorX->setValue(2, rho * u * v);
	fluxVectorX->setValue(3, rho * u * w);
	fluxVectorX->setValue(4, u * (p+e));

	// assemble flux vector in y direction:
	fluxVectorY->setValue(0, rho * v);
	fluxVectorY->setValue(1, rho * v*u);
	fluxVectorY->setValue(2, rho * v*v + p);
	fluxVectorY->setValue(3, rho * v * w);
	fluxVectorY->setValue(4, v * (p+e));

	// assemble flux vector in z direction:
	fluxVectorZ->setValue(0, rho * w);
	fluxVectorZ->setValue(1, rho * w*u);
	fluxVectorZ->setValue(2, rho * w*v);
	fluxVectorZ->setValue(3, rho * w*w + p);
	fluxVectorZ->setValue(4, w * (p+e));
};

void getEulerFluxFromPrimitiveVariables(TensorO1<double> *primitiveVariableVector, TensorO1<double> *fluxVectorX, TensorO1<double> *fluxVectorY, TensorO1<double> *fluxVectorZ ){

	TensorO1<double> consVarVect(5);
	getConservedVariableVector(primitiveVariableVector, &consVarVect);
	getEulerFluxFromConservedVariables(&consVarVect, fluxVectorX, fluxVectorY, fluxVectorZ);

};

void findRoeAverages(TensorO1<double> *leftVector, TensorO1<double> *rightVector, TensorO1<double> *avgVector){
	// Input: left and right state vectors of conserved variables
	// Output: Roe average vector of conserved variables
	// Assuming both left and right vectors are made up of conserved variables;

	double rhol, ul, vl, wl, pl, el, hl; //left state
	double rhor, ur, vr, wr, pr, er, hr; //right state
	double rhoa, ua, va, wa, pa, ea, ha; //average state
	double GAMMA = 1.4;


	//get primitive variables from the conserved variables
	TensorO1<double> leftPrimitiveVector(5);
	getPrimitiveVariableVector(leftVector, &leftPrimitiveVector);

	TensorO1<double> rightPrimitiveVector(5);
	getPrimitiveVariableVector(rightVector, &rightPrimitiveVector);

	//extract values:
	rhol = leftPrimitiveVector.getValue(0);     	rhor  = rightPrimitiveVector.getValue(0);
	ul   = leftPrimitiveVector.getValue(1);          ur   = rightPrimitiveVector.getValue(1);
	vl   = leftPrimitiveVector.getValue(2);          vr   = rightPrimitiveVector.getValue(2);
	wl   = leftPrimitiveVector.getValue(3);          wr   = rightPrimitiveVector.getValue(3);
	pl   = leftPrimitiveVector.getValue(4);          pr   = rightPrimitiveVector.getValue(4);

	//extract energy values from conserved vector
	el   = leftVector->getValue(4);			er   = rightVector->getValue(4);	//internal energy
	//find specific enthalpy
	hl   = (el + pl)/rhol;			hr   = (er + pr)/rhor;	//specific total enthalpy

	// Roe averaging
	// Refer I do like CFD (Katate Masatsuka), vol.1, page number 92, equations 3.6.112 - 3.6.115
	rhoa = sqrt(rhol * rhor);						//density
	ua = (sqrt(rhol)*ul + sqrt(rhor)*ur)/(sqrt(rhol) + sqrt(rhor));		//u velocity
	va = (sqrt(rhol)*vl + sqrt(rhor)*vr)/(sqrt(rhol) + sqrt(rhor));		//v velocity
	wa = (sqrt(rhol)*wl + sqrt(rhor)*wr)/(sqrt(rhol) + sqrt(rhor));		//w velocity
	ha = (sqrt(rhol)*hl + sqrt(rhor)*hr)/(sqrt(rhol) + sqrt(rhor));		//specific total enthalpy
	// find total internal energy for the Roe average state
	ea = 1.0/GAMMA * (ha * rhoa  + (GAMMA - 1)/2.0 * rhoa * (ua*ua + va*va + wa*wa));

	//assemble the Roe average Vector (remember this will be conserved variable vector)
	avgVector->setValue(0, rhoa);
	avgVector->setValue(1, rhoa*ua);
	avgVector->setValue(2, rhoa*va);
	avgVector->setValue(3, rhoa*wa);
	avgVector->setValue(4, ea);
};

void findEigenValuesForEuler(TensorO1<double> *conservedVariableVector, TensorO1<double> *normal, TensorO1<double> *Lambda){
	//For eigenstructure of 3D Euler's equations, refer I do like CFD, vol.1, page 77, section 3.6

	// Input: conservedVariableVector[5], normal[3];
	// output: Lambda[5];


	double rho, p;   //density and pressure
	double u, v, w;  //x,y,z components of velocity vector
	double c;	 //acoustic speed 

	// extract values of primitive variables:
	TensorO1<double> primitiveVariableVector(5);
	getPrimitiveVariableVector(conservedVariableVector, &primitiveVariableVector);

	rho = primitiveVariableVector.getValue(0) + NANO;
	u   = primitiveVariableVector.getValue(1);
	v   = primitiveVariableVector.getValue(2);
	w   = primitiveVariableVector.getValue(3);
	p   = primitiveVariableVector.getValue(4);

	// get speed of sound:
	c   = sqrt(GAMMA* p/rho);

	double un; 	  //velocity component along the given normal direction
	TensorO1<double> vel(3);
	vel.setValue(0, u);
	vel.setValue(1, v);
	vel.setValue(2, w);

	un = Math::dot(normal, &vel);

	Lambda->setValue(0, un - c);
	Lambda->setValue(1, un);
	Lambda->setValue(2, un + c);
	Lambda->setValue(3, un);
	Lambda->setValue(4, un);
};

void findEigenVectorsForEuler(TensorO1<double> *conservedVariableVector, TensorO1<double> *normal, TensorO1<double> *tangent1, TensorO1<double> *tangent2, Matrix<double> *Rn, Matrix<double> *Ln){
	//For eigenstructure of 3D Euler's equations, refer [1]:I do like CFD, vol.1, page 77, section 3.6

	// Input: conservedVariableVector[5], normal[3];
	// output: matrix Rn[5][5] and Ln[5][]

	double K = GAMMA-1.0;   // as per the reference [1]

	double rho, p;   //density and pressure
	double u, v, w;  //x,y,z components of velocity vector
	double c;	 //acoustic speed 
	double e, H;	 //total internal energy and specific enthalpy

	// extract values of primitive variables:
	TensorO1<double> primitiveVariableVector(5);
	getPrimitiveVariableVector(conservedVariableVector, &primitiveVariableVector);

	rho = primitiveVariableVector.getValue(0) + NANO;
	u   = primitiveVariableVector.getValue(1);
	v   = primitiveVariableVector.getValue(2);
	w   = primitiveVariableVector.getValue(3);
	p   = primitiveVariableVector.getValue(4);

	e   = conservedVariableVector->getValue(4);

	H   = (e + p )/ rho;

	// get speed of sound:
	c   = sqrt(GAMMA *p/rho);


	// separate components of tangents and normals
	double nx,ny,nz;
	double lx,ly,lz;
	double mx,my,mz;

	TensorO1<double> vel(3);

	vel.setValue(0,u);
	vel.setValue(1,v);
	vel.setValue(2,w);

	nx =   normal->getValue(0);  ny =   normal->getValue(1);   nz =   normal->getValue(2);
	lx = tangent1->getValue(0);  ly = tangent1->getValue(1);   lz = tangent1->getValue(2);
	mx = tangent2->getValue(0);  my = tangent2->getValue(1);   mz = tangent2->getValue(2);

	// compute normal and tangent velocity components

	double qn; 	  //velocity component along the given normal direction
	qn = Math::dot(normal, &vel);

	double ql;	  //velocity component along direction of tangent 1
	ql = Math::dot(tangent1, &vel);

	double qm; 	  //velocity component along direction of tangent 2
	qm = Math::dot(tangent2, &vel);

	double q2;	  // q_squared
	q2 = u*u + v*v + w*w;


	double mat[5][5]; //temporary matrix for right and left eigenvectors

	// compute the matrix of the right eigenvectors
	mat[0][0] =  1.0    ;   mat[0][1] = 1.0  ;   mat[0][2] =     1.0 ;   mat[0][3] = 0.0 ;   mat[0][4] = 0.0 ;   
	mat[1][0] = u - c*nx;   mat[1][1] =  u   ;   mat[1][2] = u + c*nx;   mat[1][3] =  lx ;   mat[1][4] =  mx ;   
	mat[2][0] = v - c*ny;   mat[2][1] =  v   ;   mat[2][2] = v + c*ny;   mat[2][3] =  ly ;   mat[2][4] =  my ;   
	mat[3][0] = w - c*nz;   mat[3][1] =  w   ;   mat[3][2] = w + c*nz;   mat[3][3] =  lz ;   mat[3][4] =  mz ;   
	mat[4][0] = H - qn*c;   mat[4][1] =q2/2.0;   mat[4][2] = H + qn*c;   mat[4][3] =  ql ;   mat[4][4] =  qm ;   


	// pushback to Rn
	for (Index i=0; i<5; i++){
		for (Index j=0; j<5; j++){
			Rn->setValue(i,j, mat[i][j]);
		};
	};


	// compute the matrix of left eigenvectors

	// first row
	mat[0][0] =  (K*q2/(4.0*c*c) + qn /(2.0 *c));
	mat[0][1] = -(K*u/ (2.0*c*c) + nx /(2.0 *c));
	mat[0][2] = -(K*v/ (2.0*c*c) + ny /(2.0 *c));
	mat[0][3] = -(K*w/ (2.0*c*c) + nz /(2.0 *c));
	mat[0][4] =            K/(2.0*c*c)          ;

	// second row
	mat[1][0] = 1 - (K*q2/(2.0*c*c));
	mat[1][1] =      K*u/(c*c)      ;
	mat[1][2] =      K*v/(c*c)      ;
	mat[1][3] =      K*w/(c*c)      ;
	mat[1][4] =      -K /(c*c)      ;

	// third row
	mat[2][0] =  (K*q2/(4.0*c*c) - qn /(2.0 *c));
	mat[2][1] = -(K*u/ (2.0*c*c) - nx /(2.0 *c));
	mat[2][2] = -(K*v/ (2.0*c*c) - ny /(2.0 *c));
	mat[2][3] = -(K*w/ (2.0*c*c) - nz /(2.0 *c));
	mat[2][4] =            K/(2.0*c*c)          ;

	// fourth row
	mat[3][0] = -ql ;
	mat[3][1] =  lx ;
	mat[3][2] =  ly ;
	mat[3][3] =  lz ;
	mat[3][4] = 0.0 ;

	// fifth row
	mat[4][0] = -qm ;
	mat[4][1] =  mx ;
	mat[4][2] =  my ;
	mat[4][3] =  mz ;
	mat[4][4] = 0.0 ;

	// pushback to Ln
	for (int i=0; i<5; i++){
		for (int j=0; j<5; j++){
			Ln->setValue(i,j, mat[i][j]);
		};
	};
};


void getAnalyticalShockTube(TensorO2<double> *SodAnalytical ){
	// This function returns the analytical value of the shock tube problem for given time 

	// Shock tube:
	/*
		              <-----  E.F     E.F --->   C.D -->    Shock ----> 
	          _____________________*_______*__________*___________*________
	               (2)                         (4)         (3)         (1)   

		       (2)
	          _____________________
		                       \
				        \
					 \
					  \
					   \
					    \
					     \
					      \____(4)____
					      		  |
							  |
					      		  |
							  ____(3)____					       
			                                             |
								     |   
							             |___(1)____



	The normal shock relations (and other gasdynamics relations) taken from:
	https://www.grc.nasa.gov/www/k-12/airplane/normal.html
	http://www.astro.uu.se/~hoefner/astro/teach/ch10.pdf

	*/

	double rho1, rho2, rho3, rho4;	     // density
	double u1,   u2,   u3,   u4;         // x velocity
	double v1,   v2,   v3,   v4;         // y velocity
	double w1,   w2,   w3,   w4;         // z velocity
	double p1,   p2,   p3,   p4;         // pressure
	double T1,   T2,   T3,   T4;         // temperature
	double a1,   a2,   a3,   a4;         // acoustic speed
	double M1,   M2,   M3,   M4;	     // Mach number

	double R = 287.0;

	// Right hand side  (1)
	p1 = 10.0;
	T1 = 300.0;
	rho1 = p1/(0.287*T1);	
	u1 = 0.0;
	v1 = 0.0;
	w1 = 0.0;

	// Left hand side (2)
	p2 = 50.0;
	T2 = 300.0;
	rho2 = p2/(0.287*T2);	
	u2 = 0.0;
	v2 = 0.0;
	w2 = 0.0;


	// Take a guess value of the intermediate pressure
	p3 = p2/2.0;
	p4 = p3;
	// Set all the velocities to zero (initial)
	u1 = 0.0; u2 = 0.0; u3 = 0.0; u4 = 0.0;
	v1 = 0.0; v2 = 0.0; v3 = 0.0; v4 = 0.0;
	w1 = 0.0; w2 = 0.0; w3 = 0.0; w4 = 0.0;

	double v_diff = 2.0;	// Error in velocity (difference across contact)
	double Ws; //shock speed

	while (v_diff >= MICRO){
		p4 = p3;
		// Consider stationary shock and moving reference frame
		// pressure ratio across shock
		double pr = p3/p1;
		M1 = sqrt( (pr + 0.4/2.4) * 2.4/(2*1.4));
		M3 = sqrt((0.4*M1*M1 +2.0) / (2.0 * 1.4*M1*M1 - 0.4) );
		// density ratio across shock
		double rhor=(2.4*M1*M1)/(2.0+(0.4)*M1*M1);
		// temperature ratio across shock
		double tr = ((2.0*1.4*M1*M1 - 0.4)*(0.4 * M1*M1 + 2.0)) / (2.4*2.4 * M1*M1);

		T3 = T1 * tr;
		rho3 = rho1 * rhor;
		
		// acoustic speeds
		a1 = sqrt(1.4*287*T1);
		a3 = sqrt(1.4*287*T3);

		// Consider moving shock and stationary reference frame
		// Shock speed
		Ws = M1 * a1;
		u3 = Ws - a3*M3;

		// Now consider the left hand side of the tube
		double prLeft = p4/p2;
		rho4 = rho2 * pow(prLeft, 1.0/1.4);
		T4 = p4/(rho4*0.287);

		a4 = sqrt(1.4*287*T4);
		a2 = sqrt(1.4*287*T2);

		u4 = 2.0/(0.4) * (a2 - a4);
		
	        M2 = u2/a2;	//zero
		M4=  u4/a4;
		
		v_diff=u3-u4;
		
		if (abs(v_diff)>0)
			{	p3=p3-0.0001;
			}		
		else
			{	p3=p3+0.0001;
			}
	};
	p4=p3;	

	// Values writen from left to right
	SodAnalytical->setValue(0, 0, rho2);
	SodAnalytical->setValue(1, 0, u2);
	SodAnalytical->setValue(2, 0, v2);
	SodAnalytical->setValue(3, 0, w2);
	SodAnalytical->setValue(4, 0, p2);
	SodAnalytical->setValue(5, 0, T2);
	SodAnalytical->setValue(6, 0, a2);
	SodAnalytical->setValue(7, 0, M2);

	SodAnalytical->setValue(0, 1, rho4);
	SodAnalytical->setValue(1, 1, u4);
	SodAnalytical->setValue(2, 1, v4);
	SodAnalytical->setValue(3, 1, w4);
	SodAnalytical->setValue(4, 1, p4);
	SodAnalytical->setValue(5, 1, T4);
	SodAnalytical->setValue(6, 1, a4);
	SodAnalytical->setValue(7, 1, M4);

	SodAnalytical->setValue(0, 2, rho3);
	SodAnalytical->setValue(1, 2, u3);
	SodAnalytical->setValue(2, 2, v3);
	SodAnalytical->setValue(3, 2, w3);
	SodAnalytical->setValue(4, 2, p3);
	SodAnalytical->setValue(5, 2, T3);
	SodAnalytical->setValue(6, 2, a3);
	SodAnalytical->setValue(7, 2, M3);

	SodAnalytical->setValue(0, 3, rho1);
	SodAnalytical->setValue(1, 3, u1);
	SodAnalytical->setValue(2, 3, v1);
	SodAnalytical->setValue(3, 3, w1);
	SodAnalytical->setValue(4, 3, p1);
	SodAnalytical->setValue(5, 3, T1);
	SodAnalytical->setValue(6, 3, a1);
	SodAnalytical->setValue(7, 3, M1);

};

double expp(double x,double y)
{	double z;
	z=1.0-(0.4/2*(x/y));
	return z;
}	

void getNSFluxFromConservedVariables(TensorO1<double> *conservedVariableVector, TensorO1<double> *fluxVectorX, TensorO1<double> *fluxVectorY, TensorO1<double> *fluxVectorZ, double nu ){

};

//**********************************************************************************************************************//
// TESTS, Tests
//**********************************************************************************************************************//

namespace Test{
	void TestGasdynamicsFunctions(){

		cout << "\nRunning tests on Euler flux and variable functions...\n";
		for (Index i=0; i<100; i++){
			// Testing getConservedVariableVector() function 
			TensorO1<double> primitiveVariableVector(5);

			double rho = getRandom(-100.0, 100.0);
			double u = getRandom(-100.0, 100.0);
			double v = getRandom(-100.0, 100.0);
			double w = getRandom(-100.0, 100.0);
			double p = getRandom(-100.0, 100.0);

			primitiveVariableVector.setValue(0, rho);
			primitiveVariableVector.setValue(1, u);
			primitiveVariableVector.setValue(2, v);
			primitiveVariableVector.setValue(3, w);
			primitiveVariableVector.setValue(4, p);

			TensorO1<double> conservedVariableVector(5);

			getConservedVariableVector(&primitiveVariableVector, &conservedVariableVector);

			double U[5];
			U[0] = rho;
			U[1] = rho*u;
			U[2] = rho*v;
			U[3] = rho*w;

			double GAMMA = 1.4;
			double e = p/(GAMMA -1 )  + 1.0/2.0 * rho * (u*u + v*v + w*w);
			U[4] = e;

			for (Index j=0; j<5; j++){
				assert(conservedVariableVector.getValue(j) == U[j]);
			};

			TensorO1<double>pvv(5);

			// Testing getPrimitiveVariableVector() function
			getPrimitiveVariableVector(&conservedVariableVector, &pvv);


			TensorO1<double> error(5);
			error.setValue(0, abs(pvv.getValue(0) - rho));
			error.setValue(1, abs(pvv.getValue(1) - u));
			error.setValue(2, abs(pvv.getValue(2) - v));
			error.setValue(3, abs(pvv.getValue(3) - w));
			error.setValue(4, abs(pvv.getValue(4) - p));

			assert(error.getL1Norm() < pow(10,-3) and error.getL2Norm() < pow(10,-3)and  error.getLinfNorm() < pow(10,-3) );


			// Testing Euler's flux function
			TensorO1<double> Fx(5);
			TensorO1<double> Fy(5);
			TensorO1<double> Fz(5);

			getEulerFluxFromConservedVariables(&conservedVariableVector, &Fx, &Fy, &Fz);

			rho = conservedVariableVector.getValue(0);
			u = conservedVariableVector.getValue(1)/(rho+NANO);
			v = conservedVariableVector.getValue(2)/(rho+NANO);
			w = conservedVariableVector.getValue(3)/(rho+NANO);
			e = conservedVariableVector.getValue(4);

			p = (GAMMA - 1) * (e  -   (1.0/2.0 * rho *(u*u + v*v + w*w) ) );

			// testing x flux
			double flux[5];
			flux[0] = rho*u;
			flux[1] = rho*u*u + p;
			flux[2] = rho*u*v;
			flux[3] = rho*u*w;
			flux[4] = u * (p+e);


			for (Index j=0; j<5; j++){
				assert(abs(Fx.getValue(j) - flux[j]) < SMALL);
			};

			// testing y flux
			flux[0] = rho*v;
			flux[1] = rho*v*u ;
			flux[2] = rho*v*v + p;
			flux[3] = rho*v*w;
			flux[4] = v * (p+e);


			for (Index j=0; j<5; j++){
				assert(abs(Fy.getValue(j) - flux[j]) < SMALL);
			};
			// testing z flux
			flux[0] = rho*w;
			flux[1] = rho*w*u ;
			flux[2] = rho*w*v;
			flux[3] = rho*w*w + p;
			flux[4] = w * (p+e);


			for (Index j=0; j<5; j++){
				assert(abs(Fz.getValue(j) - flux[j]) < SMALL);
			};

			// Testing primitive variable vector
			getPrimitiveVariableVector(&conservedVariableVector, &primitiveVariableVector);

			rho = primitiveVariableVector.getValue(0);
			u = primitiveVariableVector.getValue(1);
			v = primitiveVariableVector.getValue(2);
			w = primitiveVariableVector.getValue(3);
			p = primitiveVariableVector.getValue(4);

			getEulerFluxFromPrimitiveVariables(&primitiveVariableVector, &Fx, &Fy, &Fz);

			flux[0] = rho*u;
			flux[1] = rho*u*u + p;
			flux[2] = rho*u*v;
			flux[3] = rho*u*w;
			flux[4] = u * (p+e);


			for (Index j=0; j<5; j++){
				assert(abs(Fx.getValue(j) - flux[j]) < SMALL);
			};

			// testing y flux
			flux[0] = rho*v;
			flux[1] = rho*v*u ;
			flux[2] = rho*v*v + p;
			flux[3] = rho*v*w;
			flux[4] = v * (p+e);


			for (Index j=0; j<5; j++){
				assert(abs(Fy.getValue(j) - flux[j]) < SMALL);
			};
			// testing z flux
			flux[0] = rho*w;
			flux[1] = rho*w*u ;
			flux[2] = rho*w*v;
			flux[3] = rho*w*w + p;
			flux[4] = w * (p+e);


			for (Index j=0; j<5; j++){
				assert(abs(Fz.getValue(j) - flux[j]) < SMALL);
			};

		};
		cout << "Tests passed.\n";

		cout << "\nRunning Tests on Roe functions...\n";

		// tests for known values:
		TensorO1<double> leftVector(5);
		TensorO1<double> leftCVector(5);
		TensorO1<double> rightVector(5);
		TensorO1<double> rightCVector(5);
		TensorO1<double> avgVector(5);
		for (Index i=0; i<100; i++){
			leftVector.setValue(0, getRandom(1.0,10.0));
			leftVector.setValue(1, getRandom(-10.0,10.0));
			leftVector.setValue(2, getRandom(-10.0,10.0));
			leftVector.setValue(3, getRandom(-10.0,10.0));
			leftVector.setValue(4, getRandom(1.0,10.0));
			
			double rhol = leftVector.getValue(0);
			double ul = leftVector.getValue(1);
			double vl = leftVector.getValue(2);
			double wl = leftVector.getValue(3);
			double pl = leftVector.getValue(4);

			getConservedVariableVector(&leftVector, &leftCVector);


			rightVector.setValue(0, getRandom(1.0,10.0));
			rightVector.setValue(1, getRandom(-10.0,10.0));
			rightVector.setValue(2, getRandom(-10.0,10.0));
			rightVector.setValue(3, getRandom(-10.0,10.0));
			rightVector.setValue(4, getRandom(1.0,10.0));

			double rhor = rightVector.getValue(0);
			double ur = rightVector.getValue(1);
			double vr = rightVector.getValue(2);
			double wr = rightVector.getValue(3);
			double pr = rightVector.getValue(4);

			getConservedVariableVector(&rightVector, &rightCVector);

			findRoeAverages(&leftCVector, &rightCVector, &avgVector);

			double rhoa = avgVector.getValue(0);
			double ua = avgVector.getValue(1)/(rhoa +NANO);
			double va = avgVector.getValue(2)/(rhoa +NANO);
			double wa = avgVector.getValue(3)/(rhoa +NANO);
			double ea = avgVector.getValue(4);
			double pa = (GAMMA -1) * (ea - 1./2.*rhoa*(ua*ua + va*va + wa*wa));

			double rhoavg = sqrt(rhol*rhor);
			assert(Math::approx(rhoavg,rhoa));

			double uavg = (sqrt(rhol)*ul + sqrt(rhor)*ur)/(sqrt(rhol)+sqrt(rhor)); 
			assert(Math::approx(uavg,ua));

			double vavg = (sqrt(rhol)*vl + sqrt(rhor)*vr)/(sqrt(rhol)+sqrt(rhor)); 
			assert(Math::approx(vavg,va));

			double wavg = (sqrt(rhol)*wl + sqrt(rhor)*wr)/(sqrt(rhol)+sqrt(rhor)); 
			assert(Math::approx(wavg,wa));

			double el = pl/(GAMMA -1) + 1.0/2.0*rhol*(ul*ul + vl*vl + wl*wl);
			double hl = (el + pl)/rhol;

			double er = pr/(GAMMA -1) + 1.0/2.0*rhor*(ur*ur + vr*vr + wr*wr);
			double hr = (er + pr)/rhor;

			double ha;
			ha = (ea*GAMMA  - (GAMMA - 1)/2.0 * rhoa * (ua*ua + va*va + wa*wa) )/rhoa;
			double havg = (sqrt(rhol)*hl + sqrt(rhor)*hr)/(sqrt(rhol)+sqrt(rhor));

			assert(Math::approx(ha,havg));
		};
		cout << "Tests passed.\n";
		cout << endl;
	};
};
