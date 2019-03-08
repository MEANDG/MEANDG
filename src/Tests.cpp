/** @file */
#include "Tests.h"

using namespace Test;
/// each class file has a unit test written at the end under the namespace Test

/// This file runs All/ partial tests on different modules
void runAllTests(){


	printLine("~");
	std::cout << "Running ALL unit tests. This may take a while... " << std::endl;
	printLine("_");




	runMathFunctionsTests();
	runDataStructureTests();
	runGeometryTests();
	runFunctionalSpaceTests();
	runGasdynamicsFunctionsTests();




	std::cout << endl;
	printLine("~");
	std::cout << "All tests passed." << std::endl;
	printLine("~");
	std::cout << endl;

};


/** This file runs tests on all the math functions i.e. the mathFunctions.h file **/
void runMathFunctionsTests(){
	printLine("_");
	std::cout << "Running unit tests on Mathematical functions described in Math namespace (mathFunctions.h)" << std::endl;
	printLine("_");

	TestDot();
};

/** This file runs tests only on the datastructure i.e. following classes:
 *1. Point
 *2. Cell
 *3. Face
 *4. Tensor datastructure
 **/
void runDataStructureTests(){
	printLine("_");
	std::cout << "Running unit tests on data structure (Point, Cell, Face and Tensor classes)." << std::endl;
	printLine("_");


	TestTensor();
	TestPoint();
	TestFace();
	TestCell();
	TestReferenceCells();
	TestReferenceCellsAndFaces();
	TestCellAndFaceRelationship();
	TestGlobalQuadPointsMapping();
};


/** This file runs test on the geometry module. **/
void runGeometryTests(){
	printLine("_");
	std::cout << "Running unit tests on GeometryIO module." << std::endl;
	printLine("_");

	TestGeometryIO();
	TestBoundaryConditions();
};

/** This file runs test on the FunctionalSpace module. **/
void runFunctionalSpaceTests(){
	printLine("_");
	std::cout << "Running unit tests on FunctionalSpace module." << std::endl;
	printLine("_");

	//TestLagrangePolynomial();
	TestAssignment1();
	TestLGLRootsAndWeights();

	TestLagrangeMatrix();

	TestCellMatrices();
};


void runGasdynamicsFunctionsTests(){
	printLine("_");
	std::cout << "Running unit tests on Gasdynamics module." << std::endl;
	printLine("_");

	TestGasdynamicsFunctions();
};






//***********************************************************************************************************************//
// TESTS, Tests
// All the remaining tests (i.e. not declared in other classes ) are described here. 
//***********************************************************************************************************************//

namespace Test{

	void TestReferenceCellsAndFaces(){
		cout <<"\nTesting cells and reference-cell connection...\n";
		string path = "app/Tests/1x1x1";
		DG application(path);
		// Define Functional Details regarding order and roots of the functions
		int order=1; 	   //order of polynomials 
		IntFlag::intflag intType = IntFlag::inexact;  //intType: 0=inexact integration, 1 = exact integration
		
		application.setFunctionalDetails(order, intType, Solver::scalarTransportSolver, System::scalarLinear);
		//application.printFunctionalDetails();
		
		// Read the geometry files and print the details
		application.createDomain(); //now includes Jacobian, inverse Jacobian  as well as refCell array computation

		for (int i=0 ; i<4; i++){
			assert(application.referenceCells[i].getCellType() == CellType::cellType(i));
		};

		assert(application.cells[0].getCellType() == application.referenceCells[0].getCellType());

		// Printing all the reference faces
		for (int i=0 ; i<2; i++){
			assert(application.referenceFaces[i].getFaceType() == FaceType::faceType(i));
		};

		cout <<"Test::TestReferenceCellsAndFaces (in Tests.cpp) passed.\n" << endl;
	};


	void TestCellAndFaceRelationship(){
		cout << "Testing Cell and Face connectivity...\n";
		string path = "app/Tests/pitzDaily";
		DG application(path);
		// Define Functional Details regarding order and roots of the functions
		int order = 1;
		IntFlag::intflag intType = IntFlag::exact;  //intType: 0=inexact integration, 1 = exact integration
		

		application.setFunctionalDetails(order,intType, Solver::scalarTransportSolver, System::scalarLinear);
		
		// Read the geometry files and print the details
		application.createDomain(); //now includes Jacobian, inverse Jacobian  as well as refCell array computation


		 for (int i=0; i<application.noOfCells; i++){
			 for (int j=0; j<application.cells[i].getNoOfFaces(); j++){
				 if (application.cells[i].faceRelationshipArray.getValue(j) == 0){
					 assert(application.cells[i].getDefiningFace(j)->getOwnerCell()->getId() == 
					        application.cells[i].getId());
				 }
				 else{
					 assert(application.cells[i].getDefiningFace(j)->getNeighbourCell()->getId() == 
					        application.cells[i].getId());
				 }
			 };
		 };
		cout <<"Test::TestCellAndFaceRelationship (in Tests.cpp) passed.\n";
	};

	void TestGlobalQuadPointsMapping(){

		cout << "\nTesting Cell and Face Quadrature points and DOF points mapping...\n";
		// This unit test tests the mapping of the standard hex cell onto the physical hex cell (at all the 3D volume quadrature points)
		// as well as mapping of the standard quadrilateral face onto the physical quad face (at all the 2D surface quadrature points)
		// It prints the global (x,y,z) location of all the quadrature points of the required cell and the face.
		
		
		// Location of the application folder
		string path = "app/Tests/pitzDaily";


		for (Index order = 1; order < 8; order ++){
			DG application(path);
			// Define Functional Details regarding order and roots of the functions
			IntFlag::intflag intType = IntFlag::exact;  //intType: 0=inexact integration, 1 = exact integration
			

			application.setFunctionalDetails(order,intType, Solver::scalarTransportSolver, System::scalarLinear);
			
			// Read the geometry files and print the details
			application.createDomain(); //now includes Jacobian, inverse Jacobian  as well as refCell array computation

			Cell *c;
			Face *f;

			int nQuadCell;
			int nQuadFace;
			double maxError = 0;

			// Testing Quadrature points

			for (Index faceNo=0; faceNo<application.noOfFaces; faceNo++){
				f = &application.faces[faceNo];

				nQuadFace = f->getNQuad();


				for (Index QP=0; QP<nQuadFace; QP++){
					double xf, yf, zf;
					xf = f->getQuadPointsGlobalLocation()->getValue(QP,0);
					yf = f->getQuadPointsGlobalLocation()->getValue(QP,1);
					zf = f->getQuadPointsGlobalLocation()->getValue(QP,2);

					double xc, yc, zc;

					// Owner Cell

					int QPc = f->getOwnerQuadPointsArray()->getValue(QP);
					xc = f->getOwnerCell()->getQuadPointsGlobalLocation()->getValue(QPc, 0);
					yc = f->getOwnerCell()->getQuadPointsGlobalLocation()->getValue(QPc, 1);
					zc = f->getOwnerCell()->getQuadPointsGlobalLocation()->getValue(QPc, 2);

					assert(abs(xc-xf) < SMALL);
					assert(abs(yc-yf) < SMALL);
					assert(abs(zc-zf) < SMALL);

					// Neighbour cell

					if (f->getNeighbourCell() != NULL){
						QPc = f->getNeighbourQuadPointsArray()->getValue(QP);
						xc = f->getNeighbourCell()->getQuadPointsGlobalLocation()->getValue(QPc, 0);
						yc = f->getNeighbourCell()->getQuadPointsGlobalLocation()->getValue(QPc, 1);
						zc = f->getNeighbourCell()->getQuadPointsGlobalLocation()->getValue(QPc, 2);
					};


					assert(abs(xc-xf) < SMALL);
					assert(abs(yc-yf) < SMALL);
					assert(abs(zc-zf) < SMALL);
				};
							
			};

			int nDOFCell;
			int nDOFFace;

			// Testing DOF points

			for (Index faceNo=0; faceNo<application.noOfFaces; faceNo++){
				f = &application.faces[faceNo];

				nDOFFace = f->getNDOF();

				for (Index DOF=0; DOF<nDOFFace; DOF++){
					double xf, yf, zf;
					xf = f->getDOFPointsGlobalLocation()->getValue(DOF,0);
					yf = f->getDOFPointsGlobalLocation()->getValue(DOF,1);
					zf = f->getDOFPointsGlobalLocation()->getValue(DOF,2);

					double xc, yc, zc;

					// Owner Cell

					int DOFc = f->getOwnerDOFPointsArray()->getValue(DOF);
					xc = f->getOwnerCell()->getDOFPointsGlobalLocation()->getValue(DOFc, 0);
					yc = f->getOwnerCell()->getDOFPointsGlobalLocation()->getValue(DOFc, 1);
					zc = f->getOwnerCell()->getDOFPointsGlobalLocation()->getValue(DOFc, 2);


					assert(abs(xc-xf) < SMALL);
					assert(abs(yc-yf) < SMALL);
					assert(abs(zc-zf) < SMALL);

					// Neighbour Cell
					if (f->getNeighbourCell() != NULL){
						DOFc = f->getNeighbourDOFPointsArray()->getValue(DOF);
						xc = f->getNeighbourCell()->getDOFPointsGlobalLocation()->getValue(DOFc, 0);
						yc = f->getNeighbourCell()->getDOFPointsGlobalLocation()->getValue(DOFc, 1);
						zc = f->getNeighbourCell()->getDOFPointsGlobalLocation()->getValue(DOFc, 2);
					};


					assert(abs(xc-xf) < SMALL);
					assert(abs(yc-yf) < SMALL);
					assert(abs(zc-zf) < SMALL);
				};
							
			};

			printProgressBar(order*1.0/7.0 * 100.0);
		};
		cout << endl;
		cout <<"Test::TestGlobalPointsMapping (int Tests.cpp) passed.\n";
	};


	void TestCellMatrices(){
		cout << "Testing Cell Matrices...\n";

		string path = "app/Tests/1x1x1";
		ofstream file;
		string filename = "tests/CellMatrices.dat";
		file.open(filename.c_str());

		file << "_______________________________________________________________________" << endl;
		file << "Inexact Integration" << endl;
		file << "_______________________________________________________________________" << endl;

		for (Index order = 1; order < 5; order ++){
			DG application(path);
			// Define Functional Details regarding order and roots of the functions
			IntFlag::intflag intType = IntFlag::inexact;  //intType: 0=inexact integration, 1 = exact integration
			
			application.setFunctionalDetails(order, intType, Solver::scalarTransportSolver, System::scalarLinear);
			//application.printFunctionalDetails();
			
			// Read the geometry files and print the details
			application.createDomain(); //now includes Jacobian, inverse Jacobian  as well as refCell array computation

			cout << "Order: " << order << "   (Inexact integration)" << endl;
			application.computeCellMatrix();


			assert(fileExists(filename));

			file << "Order: " << order << endl;
			for (Index c=0; c<application.noOfCells; c++){
				application.cells[c].printMassMatrix(file);
				application.cells[c].printInverseMassMatrix(file);
				application.cells[c].printVandermondeMatrix(file);
				application.cells[c].printDifferentialMatrix(file);
				application.cells[c].printDTildeMatrix(file);
				application.cells[c].printFluxMatrix(file);
				application.cells[c].printFTildeMatrix(file);
			};

			file << "_______________________________________________________________________" << endl;

		};

		file << "_______________________________________________________________________" << endl;
		file << "Exact Integration" << endl;
		file << "_______________________________________________________________________" << endl;

		for (Index order = 1; order < 5; order ++){
			DG application(path);
			// Define Functional Details regarding order and roots of the functions
			IntFlag::intflag intType = IntFlag::exact;  //intType: 0=inexact integration, 1 = exact integration
			
			application.setFunctionalDetails(order, intType, Solver::scalarTransportSolver, System::scalarLinear);
			//application.printFunctionalDetails();
			
			// Read the geometry files and print the details
			application.createDomain(); //now includes Jacobian, inverse Jacobian  as well as refCell array computation

			cout << "Order: " << order << "   (Exact integration)" << endl;
			application.computeCellMatrix();


			assert(fileExists(filename));

			file << "Order: " << order << endl;
			for (Index c=0; c<application.noOfCells; c++){
				application.cells[c].printMassMatrix(file);
				application.cells[c].printInverseMassMatrix(file);
				application.cells[c].printVandermondeMatrix(file);
				application.cells[c].printDifferentialMatrix(file);
				application.cells[c].printFluxMatrix(file);
				application.cells[c].printFTildeMatrix(file);
			};

			file << "_______________________________________________________________________" << endl;

		};

		file.close();
		cout << "Test::TestCellMatrices() (in Tests.cpp) passed.\n";
		cout << "Cell matrices (for Exact and Inexact integration) are printed in tests/CellMatrices.dat. Kindly inspect.\n" << endl;
	};


	void TestBoundaryConditions(){

		cout << "Testing Boundary Conditions...\n";
		string path = "app/Tests/pitzDaily";
		DG application(path);
		// Define Functional Details regarding order and roots of the functions
		int order = 1;
		IntFlag::intflag intType = IntFlag::exact;  //intType: 0=inexact integration, 1 = exact integration
		

		application.setFunctionalDetails(order,intType, Solver::scalarTransportSolver, System::scalarLinear);
		
		// Read the geometry files and print the details
		application.createDomain(); //now includes Jacobian, inverse Jacobian  as well as refCell array computation



		
		// Define temporal and integrator details
		double deltaT = 1e-3,startTime=0.0,endTime=5,printTime=1e-1;int integratorType=3;
		application.setTemporalDetails(deltaT,startTime,endTime,printTime,integratorType);

		// Define and read the variables from 0 folder
		string VariableName[2] = {"T", "U"};
		int noOfVariable = 2;
		int VariableSize[2] = {1, 3};
		application.generateVariableArray(noOfVariable,VariableName,VariableSize);
		application.readData(0);

		assert(application.getTotNoOfVariables() == 4);

		assert(application.bcond[0].getName() == "inlet" and application.bcond[0].getStartFace() == 24170 and application.bcond[0].getNoOfFaces() == 30);
		assert(application.bcond[1].getName() == "outlet" and application.bcond[1].getStartFace() == 24200 and application.bcond[1].getNoOfFaces() == 57);
		assert(application.bcond[2].getName() == "upperWall" and application.bcond[2].getStartFace() == 24257 and application.bcond[2].getNoOfFaces() == 223);
		assert(application.bcond[3].getName() == "lowerWall" and application.bcond[3].getStartFace() == 24480 and application.bcond[3].getNoOfFaces() == 250);
		assert(application.bcond[4].getName() == "frontAndBack" and application.bcond[4].getStartFace() == 24730 and application.bcond[4].getNoOfFaces() == 24450);


		cout << "Test::TestBoundaryConditions (in Tests.cpp) passed.\n";


	};
	

	void TestScalarTransportSolver(){
		string path;
		path = "app/circAdv_8x8x8";
		DG application(path);

		// Define Functional Details regarding order and roots of the functions
		int order=1; 
		IntFlag::intflag intType = IntFlag::exact;
		
		Solver::solver solverType = Solver::scalarTransportSolver;      // refer include/Solvers.h for other solvers
		System::system systemType = System::scalarLinear;       	// refer include/Systems.h for available systems of equations
		application.setFunctionalDetails(order,intType, solverType, systemType);

		application.printFunctionalDetails();
		
		// Read the geometry files and print the details
		application.createDomain();
		application.printDomain();
		
		// Define temporal and integrator details
		double deltaT = 1e-2,startTime=0,endTime=2*PI,printTime=0.1;int integratorType=3;
		//double deltaT = 1e-4,startTime=0.0,endTime=10e-4,printTime=1;int integratorType=3;
		application.setTemporalDetails(deltaT,startTime,endTime,printTime,integratorType);
		application.printTemporalDetails();

		// Define and read the variables from 0 folder
		string VariableName[2] = {"T", "U"};
		int noOfVariable = 2;
		int VariableSize[2] = {1,3};

		application.generateVariableArray(noOfVariable,VariableName,VariableSize);
		//application.applyIC(); //reads array from 0 folder then creates a 1 folder with correct IC. copy this manually to 0 folder and change the approriate field inside of all files.
		application.readData(0); //reference solution is stored same as IC


		cout << "\n\nPrinting Boundary conditions\n";
		cout << "No. \t Name \t         Starting face \t        No. of Faces \n";
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		for (int i=0; i<application.noOfBoundaryConditions; i++){
			cout << application.bcond[i].getId() << "\t" << application.bcond[i].getName() << " \t Starting Face: " <<  application.bcond[i].getStartFace() << "\tNumber of faces: " << application.bcond[i].getNoOfFaces() << endl;
		};
		
		// Running Simulation
		application.runApplication();

		application.integrateCellVariable();
		application.writeVariableArray(endTime);

	};



	void TestEulerEquations1(){
		string path;
		path = "app/linAdvEuler";
		DG application(path);

		// Define Functional Details regarding order and roots of the functions
		int order=1; 
		IntFlag::intflag intType = IntFlag::inexact;
		
		Solver::solver solverType = Solver::LLF;      // refer include/Solvers.h for other solvers
		System::system systemType = System::Euler;       	// refer include/Systems.h for available systems of equations
		application.setFunctionalDetails(order,intType, solverType, systemType);

		application.printFunctionalDetails();
		
		// Read the geometry files and print the details
		application.createDomain();
		application.printDomain();
		
		// Define temporal and integrator details
		//double deltaT = 1e-2,startTime=0,endTime=12,printTime=0.2;int integratorType=3;
		double deltaT = 5e-5,startTime=0.0,endTime=100*deltaT,printTime=10*deltaT;int integratorType=3;
		application.setTemporalDetails(deltaT,startTime,endTime,printTime,integratorType);
		application.printTemporalDetails();

		// Define and read the variables from 0 folder
		string VariableName[3] = {"rho", "U", "p"};
		int noOfVariable = 3;
		int VariableSize[3] = {1,3,1};

		application.generateVariableArray(noOfVariable,VariableName,VariableSize);
		//application.applyIC(); //reads array from 0 folder then creates a 1 folder with correct IC. copy this manually to 0 folder and change the approriate field inside of all files.
		application.readData(startTime); //reference solution is stored same as IC
		application.writeVariableRefArray(startTime);


		cout << "\n\nPrinting Boundary conditions\n";
		cout << "No. \t Name \t         Starting face \t        No. of Faces \n";
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		for (int i=0; i<application.noOfBoundaryConditions; i++){
			cout << application.bcond[i].getId() << "\t" << application.bcond[i].getName() << " \t Starting Face: " <<  application.bcond[i].getStartFace() << "\tNumber of faces: " << application.bcond[i].getNoOfFaces() << endl;
		};

		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		cout << "Boundary type for all the variables:  \n";
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		for (int i=0; i<application.noOfBoundaryConditions; i++){
			// At a particular boundary condition
			cout << application.bcond[i].getName() << endl;
			cout << " [ ";
			for (Index v=0; v<application.getTotNoOfVariables(); v++){
				cout << application.bcond[i].getVariableTypeArray()->getValue(v) << "  ";
			};
			cout << " ]\n" << endl;
		};
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		cout << "Boundary values for all the variables:  \n";
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		for (int i=0; i<application.noOfBoundaryConditions; i++){
			cout << application.bcond[i].getName() << endl;
			cout << " [ ";
			for (Index v=0; v<application.getTotNoOfVariables(); v++){
				cout << application.bcond[i].getVariableValueArray()->getValue(v) << "  ";
			};
			cout << " ]\n" << endl;
		};
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		cout << "Boundary functions for all the variables:  \n";
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		for (int i=0; i<application.noOfBoundaryConditions; i++){
			cout << application.bcond[i].getName() << endl;
			cout << " [ ";
			for (Index v=0; v<application.getTotNoOfVariables(); v++){
				cout << application.bcond[i].getFunc2bcondTypeArray()->getValue(v) << "  ";
			};
			cout << " ]\n" << endl;

		};
		
		// Running Simulation
		application.runApplication();

		application.integrateCellVariable();
		application.writeVariableArray(endTime);
		application.writeVariableRefArray(endTime);

	};



}; //EnD NamespacE TesT
