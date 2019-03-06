#include"RiemannSolver.h"

using namespace std;

void scalarRiemannSolver(Face *face, int rkstep){
	TensorO1<double> vel(3);
	double T;

	// find the velocity normal to the face
	double un = 0.0;

	int J = 1; //some integer
	// velocity vector is not getting updated. so rkstep = 0 for this vector;
	vel.setValue(0, face->getOwnerCell()->getVariable(0, 1, face->getOwnerDOFPointsArray()->getValue(J)));
        vel.setValue(1, face->getOwnerCell()->getVariable(0, 2, face->getOwnerDOFPointsArray()->getValue(J)));
        vel.setValue(2, face->getOwnerCell()->getVariable(0, 3, face->getOwnerDOFPointsArray()->getValue(J)));


	un = Math::dot(&vel, face->getNormal());


	TensorO1<double> fluxx(face->getNDOF());
	TensorO1<double> fluxy(face->getNDOF());
	TensorO1<double> fluxz(face->getNDOF());


	if (un > 0){ 
		//fluid going out of the element through face
		for (Index DOF=0; DOF< face->getNDOF(); DOF++){
			T     = face->getOwnerCell()->getVariable(   rkstep, 0, face->getOwnerDOFPointsArray()->getValue(DOF));
			vel.setValue(0, face->getOwnerCell()->getVariable(0, 1, face->getOwnerDOFPointsArray()->getValue(DOF)));
			vel.setValue(1, face->getOwnerCell()->getVariable(0, 2, face->getOwnerDOFPointsArray()->getValue(DOF)));
			vel.setValue(2, face->getOwnerCell()->getVariable(0, 3, face->getOwnerDOFPointsArray()->getValue(DOF)));
			// Note: velocity vector is not getting updated. so rkstep = 0 for this vector;

			fluxx.setValue(DOF, vel.getValue(0) * T);
			fluxy.setValue(DOF, vel.getValue(1) * T);
			fluxz.setValue(DOF, vel.getValue(2) * T);
			// 0 indicates variable number. for scalar, it is 0.
		};
	}
	else if (un < 0){
		//fluid entering the cell
		if (face->getNeighbourCell() != NULL){
			for (int DOF=0; DOF< face->getNDOF(); DOF++){
				T     = face->getNeighbourCell()->getVariable(   rkstep, 0, face->getNeighbourDOFPointsArray()->getValue(DOF));
				vel.setValue(0, face->getNeighbourCell()->getVariable(0, 1, face->getNeighbourDOFPointsArray()->getValue(DOF)));
				vel.setValue(1, face->getNeighbourCell()->getVariable(0, 2, face->getNeighbourDOFPointsArray()->getValue(DOF)));
				vel.setValue(2, face->getNeighbourCell()->getVariable(0, 3, face->getNeighbourDOFPointsArray()->getValue(DOF)));

				fluxx.setValue(DOF, vel.getValue(0) * T);
				fluxy.setValue(DOF, vel.getValue(1) * T);
				fluxz.setValue(DOF, vel.getValue(2) * T);
			};
		}
		else{
			// Boundary conditions will overwrite these values
			for (int DOF=0; DOF< face->getNDOF(); DOF++){
				fluxx.setValue(DOF, 0.0);
				fluxy.setValue(DOF, 0.0);
				fluxz.setValue(DOF, 0.0);
			};
		}
	}
	else{   //un = 0.0
		for (int DOF=0; DOF< face->getNDOF(); DOF++){
			// no flux
			fluxx.setValue(DOF, 0.0);
			fluxy.setValue(DOF, 0.0);
			fluxz.setValue(DOF, 0.0);
		};
	};



	TensorO1<double> FLUX(3);
	// Normal flux at the face
	for (int DOF=0; DOF< face->getNDOF(); DOF++){

		// create a local flux vector at DOF point
		FLUX.setValue(0, fluxx.getValue(DOF));
		FLUX.setValue(1, fluxy.getValue(DOF));
		FLUX.setValue(2, fluxz.getValue(DOF));

		face->setFlux(0, DOF, Math::dot(&FLUX, face->getNormal()));
		// 0 in face.Flux[0][DOF]  stands for scalar variable
	};
};




void BurgerRiemannSolver(Face *face, int rkstep){
};

void centralDifferenceHeatEquation(Face *face, int rkstep){
};

void RoeSolver(Face *face, int rkstep){
};

void RoeSolverNS(Face *face, int rkstep){
};

void LLFSolver(Face *face, int rkstep){
	// Local Lax Frierich

	// Create local vectors
	// remeber: these are conserved variable vectors not primitive
	TensorO1<double> leftStateVector(5);
	TensorO1<double> rightStateVector(5);
	TensorO1<double> averageStateVector(5);
	TensorO1<double> deltaQ(5);    		//rightStateVector - leftStateVector;
	TensorO1<double> diffusion(5); 		//diffusion vector which is added to the central difference

	// Next, local flux vectors at a particular DOF location
	TensorO1<double> fluxVectorX(5);
	TensorO1<double> fluxVectorY(5);
	TensorO1<double> fluxVectorZ(5);

	TensorO1<double> fluxVectorLeft(5); 	//Normal flux based on the left variable vector
	TensorO1<double> fluxVectorRight(5); 	//and based on the right variable vector

	TensorO1<double> LambdaL(5);	     	//eigenvalues
	TensorO1<double> LambdaR(5);	     	//eigenvalues
	TensorO1<double> LambdaA(5);	     	//eigenvalues
	double LambdaMax;		//max eigenvalue
	double alpha = 0.0;

	face->getFluxArray()->setAll(0.0);

	TensorO1<double> momentum(3);
	TensorO1<double> momentumR(3);
	TensorO1<double> LSR(5);
	TensorO1<double> RSR(5);


	// Find flux vectors for LLF solver
	if (face->getNeighbourCell() != NULL){
		for (int DOF=0; DOF< face->getNDOF(); DOF++){
			int DOFc = face->getOwnerDOFPoint(DOF);
			int DOFn = face->getNeighbourDOFPoint(DOF);

			// copy values of the conserved variable vector
			for (int nVar=0; nVar<5; nVar++){
				leftStateVector.setValue(nVar, face->getOwnerCell()->getVariable(rkstep, nVar, DOFc));
				rightStateVector.setValue(nVar, face->getNeighbourCell()->getVariable(rkstep,nVar,DOFn));
			};

			// first left state
			for (int k=0; k<3; k++){
				momentum.setValue(k, leftStateVector.getValue(k+1));
			};

			face->rotateNormal(&momentum, &momentumR);

			LSR.setValue(0, leftStateVector.getValue(0));
			LSR.setValue(1, momentumR.getValue(0));
			LSR.setValue(2, momentumR.getValue(1));
			LSR.setValue(3, momentumR.getValue(2));
			LSR.setValue(4, leftStateVector.getValue(4));
			
			// right state
			for (int k=0; k<3; k++){
				momentum.setValue(k, rightStateVector.getValue(k+1));
			};

			face->rotateNormal(&momentum, &momentumR);

			RSR.setValue(0, rightStateVector.getValue(0));
			RSR.setValue(1, momentumR.getValue(0));
			RSR.setValue(2, momentumR.getValue(1));
			RSR.setValue(3, momentumR.getValue(2));
			RSR.setValue(4, rightStateVector.getValue(4));
			
			// compute eigenvalues and eigenvectors based on this state
			findEigenValuesForEuler(&leftStateVector   , face->getNormal(), &LambdaL);
			findEigenValuesForEuler(&rightStateVector  , face->getNormal(), &LambdaR);

			// filling diagAbsLambda matrix
			alpha = max(LambdaL.getLinfNorm(),  LambdaR.getLinfNorm());


			// find flux vector based on left and right states
			// refer src/gasdynamics for implementation

			// Left:
			getEulerFluxFromConservedVariables(&LSR, &fluxVectorX, &fluxVectorY, &fluxVectorZ);
			// find normal flux (only x directional flux used because the momentum vector is rotated normal to face)
			for (int k=0; k<5; k++){
				fluxVectorLeft.setValue(k, fluxVectorX.getValue(k));
			};

			// Right:
			getEulerFluxFromConservedVariables(&RSR, &fluxVectorX, &fluxVectorY, &fluxVectorZ);
			// find normal flux 
			for (int k=0; k<5; k++){
				fluxVectorRight.setValue(k, fluxVectorX.getValue(k));
			};


			// last required element is DeltaQ vector
			for (int k=0; k<5; k++){
				deltaQ.setValue(k, RSR.getValue(k) - LSR.getValue(k));
			};

			// LLF flux formula
			for (int k=0; k<5; k++){
				double FLUX = 1.0/2.0*(fluxVectorRight.getValue(k) + fluxVectorLeft.getValue(k)) - 
					      1.0/2.0 * abs(alpha)* deltaQ.getValue(k);
				face->setFlux(k,DOF,FLUX); 
			};

			// first left state
			for (int k=0; k<3; k++){
				momentumR.setValue(k, face->getFlux(k+1,DOF));
			};
			face->rotateBackNormal(&momentumR, &momentum);

			for (int k=0; k<3; k++){
				face->setFlux(k+1,DOF, momentum.getValue(k));
			};

		};

	} 	// if neighbour != NULL condition finished

	else{
		//implement Boundary conditions
		for (int DOF=0; DOF< face->getNDOF(); DOF++){
			for (int k=0; k<5; k++){
				face->setFlux(k,DOF,0.0); 
			};
		};
	};
};

void LLFSolverNS(Face *face, int rkstep){
};

