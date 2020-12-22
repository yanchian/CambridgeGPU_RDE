/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
//#define GHOST
#define REACTIVE
#pragma once
#include "core.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <deque>
#include <queue>
#include <ctime>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <signal.h>
#include "boost/thread.hpp"
#include "Matrix.hpp"
#include "MatrixOperations.hpp"

#define USE_GL

#define MATERIALS 2
#define DIMENSIONS 2
#define GC 2

#include "grid.hpp"

template<Processor P>
struct Mesh {
  typedef Mesh2D<real, P, NUMBER_VARIABLES> type;
};

template<Processor P>
struct LevelSet {
  typedef Mesh2D<real, P, 2> type;
};

StridedCell<real, NUMBER_VARIABLES> typedef Cell;

/*///////////////////////////////////////////////////////////////////////////////////////////
// Read input flow field variables
int read() {
double var[1];
std::ifstream inputfile;
inputfile.open("INPUT.txt");

while (! inputfile.eof()){
    inputfile >> var[0];
}
inputfile.close();
return 0;
};
///////////////////////////////////////////////////////////////////////////////////////////*/


const real X1 = 0.0;
const real X2 = 100.0;
const real L1 = 48.0;//192
const real L2 = 150.0;
const real PCJ = 17.75;
const real rho_CJ = 1.715;
const real pShock = 2.0 * PCJ;
const real rhoShock = 1.0 * rho_CJ;
const real pAir = 1.0;
const real rhoAir = 1.0;
const real Tol = 1.0e-6;

const real Pst = 1.0*20;
const real Pcr =  1.0*10.843;

//const real Pcr1 = 18.79;
//const real Pcr2 = 10.514;
//const real Pcr2 = 10.51;
//const real Ma1 = 2.139;
//const real Ma2 = 0.308;
const real Mcj = 5.517;
const real Tst =4.3;

// Sinusoidal perturbation
const real PI = 3.14159265359;
const real AMP = 2.5;           // Amplitude
const real LAMBDA = 120.0;       // Wavelength
const real wavelength = 200;
const real cellData = 10;//10
const real xinput = 190;
const real yinput = 192;
const real nozzleD = 10;
const real wallD= 0;
const real slotNumber = 40;
/*
// Put in a 90 degree opening
const real XCORNER = 1250.0;
const real YCORNER = 240.0;
const real GA = 0.04;
const real SPACING = 40.0;
*/

/*
// Put in a rectnagular block as a small perturbation
const real block_x1 = 1400.0;
const real block_x2 = 1450.0;
const real block_y1 = 280.0;
const real block_y2 = 285.0;
*/

/*
// Parameters for 2-step kinetics (unstable H2/O2)
const real Q = 21.365;
const real TS = 5.0373;
const real EI = 5.414 * TS;
const real ER = 1.0 * TS;
const real KI = 1.0022;
const real KR = 4.0;
const real specific_heat_ratio = 1.32;
*/

/*
// Parameters for 2-step kinetics (stable C2H2/O2/Ar)
const real Q = 19.7;
const real TS = 7.6051;
const real EI = 4.8 * TS;
const real ER = 1.0 * TS;
const real KI = 0.139;
const real KR = 0.2;
const real specific_heat_ratio = 1.212;
*/
const real Q = 25.31;
const real TS = 5.7353;
const real EI = 6.52 * TS;
const real ER = 1.0 * TS;
const real KI = 1.0538;
const real KR = 1.0; //3.74
const real specific_heat_ratio = 1.32;
// Parameters for single-step Arrhenius kinetics
//const real Q = 50;
//const real KR = 85.5; 
//const real Ea = 30.0;
//const real KR = 16.45;
//const real KR = 7.25;
//const real KR = 3.55;
//const real Ea = 20.0;
//const real KR = 8.5;
//const real KR = 1.6;
//const real Ea = 15.0;
//const real KR = 80.2;
//const real Ea = 30.0;
//const real specific_heat_ratio = 1.2;

real Counter = 1.0;
real Check = 1.0;
//const int len_inert = 500; //in terms of cells (multiplied by resolution)
real density_new = 1.0;
const real Delta_rho = 0.0;
const real MFR = 0.0;
const real frame_interface = 0.0;
const real Plot_step = 49.5;
real Plot_counter = 10.0;
const int Skip_lines = 5.0;

__device__ __host__ __forceinline__ real gamma(const Cell u) {
  return specific_heat_ratio;
  
}
__device__ __host__ __forceinline__ real gamma() {
  return specific_heat_ratio;
}

__device__ __host__ __forceinline__ real p0(const Cell u) {
  return 0.0;
  //return 0.87e8;
}
__device__ __host__ __forceinline__ real p0() {
  return 0.0;
  //return 0.87e8;
}

#include "boundaryconditions.hpp"
#include "flux.hpp"
#include "wavespeed.hpp"
#include "HLLC.hpp"
#include "Solver.hpp"
#include "initialconditions.hpp"
#include "render.hpp"
#include "opengl.hpp"
#ifdef GHOST
#include "ghostfluid.hpp"
#endif
#ifdef REACTIVE
#include "source.hpp"
#include "shockdetect.hpp"
#endif

struct ImageOutputs {
  std::string prefix;
  int plotVariable;
  ColourMode colourMode;
  real min;
  real max;
};

#include "SDF/BoundaryMesh.hpp"
#include "SDF/Polyhedron.cu"
#include "SDF/ConnectedEdge.cu"
#include "SDF/Edge.cu"
#include "SDF/Face.cu"
#include "SDF/Vertex.cu"
#include "SDF/ConnectedFace.cu"
#include "SDF/ConnectedVertex.cu"
#include "SDF/ScanConvertiblePolygon.cu"
#include "SDF/ScanConvertiblePolyhedron.cu"
#include "SDF/BoundaryMesh.cu"

#include "kernels/boundaryconditions.ipp"
#include "kernels/flux.ipp"
#ifdef GHOST
#include "kernels/ghostfluid.ipp"
#endif
#ifdef REACTIVE
#include "kernels/source.ipp"
#include "kernels/shockdetect.ipp"
#endif
#include "kernels/HLLC.ipp"
#include "kernels/wavespeed.ipp"
