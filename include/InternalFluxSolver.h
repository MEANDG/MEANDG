#ifndef INTERNALFLUXSOLVER_SOLVER
#define INTERNALFLUXSOLVER_SOLVER

#include"sysInclude.h"
#include"Face.h"
#include"Cell.h"
#include "Gasdynamics.h"
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Scalar Riemann Solvers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void scalarTransportFlux(Cell *cell, int rkstep);
void BurgerFlux(Cell *cell, int rkstep);
void HeatEquationFlux(Cell *cell, int rkstep);
void EulerFlux(Cell *cell, int rkstep);
void NavierStokesFlux(Cell *cell, int rkstep);

#endif
