/** @file */
#ifndef TESTS
#define TESTS

#include "sysInclude.h"
#include "mathFunctions.h"
#include "Point.h"
#include "Face.h"
#include "Cell.h"
#include "GeometryIO.h"
#include "FunctionalSpaces.h"
#include "RefCell.h"
#include "RefFace.h"
#include "Gasdynamics.h"
#include "DG.h"

using namespace Test;

void runAllTests();
void runDataStructureTests();
void runMathFunctionsTests();
void runGeometryTests();
void runFunctionalSpaceTests();
void runGasdynamicsFunctionsTests();



//***********************************************************************************************************************//
// TESTS, Tests
// All the remaining tests (i.e. not declared in other classes ) are described here. 
//***********************************************************************************************************************//

namespace Test{
	void TestReferenceCellsAndFaces();
	void TestGlobalQuadPointsMapping();
	void TestCellAndFaceRelationship();
	void TestCellMatrices();
	void TestBoundaryConditions();
	void TestScalarTransportSolver();
	void TestEulerEquations1();
};

#endif
