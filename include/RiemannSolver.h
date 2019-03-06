#ifndef RIEMANN_SOLVER
#define RIEMANN_SOLVER

#include "sysInclude.h"
#include "mathFunctions.h"
#include "Face.h"
#include "Cell.h"
#include "Gasdynamics.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Scalar Riemann Solvers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void scalarRiemannSolver(Face *face, int rkstep);
void BurgerRiemannSolver(Face *face, int rkstep);
void centralDifferenceHeatEquation(Face *face, int rkstep);
void RoeSolver(Face *face, int rkstep);
void RoeSolverNS(Face *face, int rkstep);
void LLFSolver(Face *face, int rkstep);
void LLFSolverNS(Face *face, int rkstep);

#endif
