/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "boundaryconditions.hpp"
#include "ghost.hpp"
#include <math.h>
template<BoundaryConditions BCs, bool XDIR, bool downstream>
__global__ void setBoundaryConditionsKernel(Mesh<GPU>::type u) {
  const bool YDIR = !XDIR, upstream = !downstream;

  const int i = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();



// for the top and bottom boundaries
  if (XDIR && u.exists(i, 0)) {
      if (downstream) {
		  //top
			for (int k = 0; k < NUMBER_VARIABLES; k++) {
				
				//u(i, -2, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 1, k);
                //u(i, -1, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 0, k);	
				
				for (int nozzlenumber = 1; nozzlenumber <= slotNumber ; nozzlenumber++) {
					
				//if (i >= (nozzlenumber * 2 - 1) * nozzleD * cellData + 1 && i <= nozzlenumber * 2 * nozzleD * cellData) {
				if (i >= ((nozzlenumber - 1) * (wallD + nozzleD) + wallD)*cellData && i <= nozzlenumber * (wallD + nozzleD) * cellData - 1) { 	
				//if (i >= 0 && i <= 3500) {	
				
				real temp[NUMBER_VARIABLES];
				real p[NUMBER_VARIABLES];
				
				const real gamma1 = specific_heat_ratio;
				
				real temperature = Tst;
				real Ma = 0;
				temp[DENSITY] = Pst/Tst;
				temp[XVELOCITY] = 0;
				temp[YVELOCITY] = 0.0;
				temp[PRESSURE] = Pst;
				temp[LAMBDA0] = 0.0;
				temp[LAMBDA1] = 0.0;
				
				//conservativeToPrimitiveInPlace(u(i, 0));
				conservativeToPrimitive(u(i, 0), p);
				
				//if (u(i, 0,PRESSURE) >= Pst) {
				if (p[PRESSURE] >= Pst) {	
				//u(i, -2, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 1, k);
                //u(i, -1, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 0, k);
				//temp[PRESSURE] = p[PRESSURE];	
				//temperature = Tst*pow((temp[PRESSURE]/Pst),(gamma1-1)/gamma1);
                //temp[DENSITY] = temp[PRESSURE]/temperature;
                //temp[DENSITY] = p[DENSITY];
						
				//temp[XVELOCITY] = 0,0;
				//temp[YVELOCITY] = 0.0;							
				//temp[LAMBDA0] = p[LAMBDA0];
			    //temp[LAMBDA1] = p[LAMBDA1];
				u(i, -1, k) = u(i, 0, k);
				u(i, -2, k) = u(i, 1, k);
				u(i, -1, XMOMENTUM) = -u(i, 0, XMOMENTUM);
				u(i, -2, XMOMENTUM) = -u(i, 1, XMOMENTUM);
			    u(i, -1, YMOMENTUM) = -u(i, 0, YMOMENTUM);
				u(i, -2, YMOMENTUM) = -u(i, 1, YMOMENTUM);

				//u(i, -1, YMOMENTUM) = 0;
				//u(i, -2, YMOMENTUM) = 0;

                //primitiveToConservativeInPlace(temp);	
				//u(i, -2, k) = temp[k] ;
				//u(i, -1, k) = temp[k] ;
				//primitiveToConservativeInPlace(temp);	
				//u(i, -2, k) = temp[k] ;
				//u(i, -1, k) = temp[k] ;
				}
				
				//else if (u(i, 0, PRESSURE) >= Pcr1 && u(i, 0, PRESSURE) < Pst) {
			    else if (p[PRESSURE] > Pcr && p[PRESSURE] < Pst) {
				//temp[YVELOCITY] = sqrt(2*(pow((Pst/u(i, 0, 3)),((gamma1-1)/gamma1))-1)/(gamma1-1));
				//temp[PRESSURE] = Pst*pow(1+(gamma1-1)/2*pow(temp[YVELOCITY],2),-gamma1/(gamma1-1));
				//temperature = Tst*pow(1+((gamma1-1)/2)*pow(temp[YVELOCITY],2),-1);
				//temp[DENSITY] = temp[PRESSURE]/temperature;
				temp[PRESSURE] = p[PRESSURE];
				temperature = Tst*pow((temp[PRESSURE]/Pst),(gamma1-1)/gamma1);
				temp[DENSITY] = temp[PRESSURE]/temperature;
				temp[XVELOCITY] = - p[XVELOCITY];
				temp[YVELOCITY] = sqrt((2*gamma1/(gamma1-1))*Tst*(1-(pow(temp[PRESSURE]/Pst,(gamma1-1)/gamma1))));
				//temp[LAMBDA0]=p[LAMBDA0];
				primitiveToConservativeInPlace(temp);	

				u(i, -2, k) = temp[k] ;
				u(i, -1, k) = temp[k] ;
				}
				
				//else if (u(i, 0, PRESSURE) >= Pcr2 && u(i, 0, PRESSURE) < Pcr1) {
					/*
				else if (p[PRESSURE] >= Pcr2 && p[PRESSURE] < Pcr1) {
				temp[PRESSURE]= p[PRESSURE];
				Ma = sqrt(-1/(gamma1-1)+sqrt(pow(1/(gamma1-1),2)+(2/(gamma1-1))*pow((2/(gamma1+1)),(gamma1+1)/(gamma1-1))*pow(Pst/temp[PRESSURE],2)*pow(0.5,2)));
				temperature = Tst*pow(1+((gamma1-1)/2)*pow(Ma,2),-1);
				temp[YVELOCITY] = sqrt(gamma1*temperature) * Ma;
				//temp[YVELOCITY] = sqrt(-1/(gamma1-1)+sqrt(pow(1/(gamma1-1),2)+(2/(gamma1-1))*pow((2/(gamma1+1)),(gamma1+1)/(gamma1-1))*pow(Pst/u(i, 0, 3),2)*pow(0.5,2)));
				temp[DENSITY] = temp[PRESSURE]/temperature;
				//temp[LAMBDA0]=p[LAMBDA0];
				primitiveToConservativeInPlace(temp);	
				u(i, -2, k) = temp[k] ;
				u(i, -1, k) = temp[k] ;
				}
				*/
				//if (u(i, 0, PRESSURE) < Pcr2) {		
				else
				{				
				//temp[YVELOCITY] = Ma1;
				temp[PRESSURE] = Pcr;
				temp[YVELOCITY] = sqrt((2*gamma1/(gamma1-1))*Tst*(1-(pow(temp[PRESSURE]/Pst,(gamma1-1)/gamma1))));				
				temperature = Tst*pow((temp[PRESSURE]/Pst),(gamma1-1)/gamma1);
				temp[DENSITY] = temp[PRESSURE]/temperature;		
				temp[XVELOCITY] = - p[XVELOCITY];
				//temp[LAMBDA0]=p[LAMBDA0];
				primitiveToConservativeInPlace(temp);	
				u(i, -2, k) = temp[k] ;
				u(i, -1, k) = temp[k] ;
				}


   //      u(i, -2, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 1, k);
   //      u(i, -1, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 0, k);
   // for periodic BCs
   //     u(i, -2, k) = u(i, u.activeNy() - 2, k); 
   //     u(i, -1, k) = u(i, u.activeNy() - 1, k);
   // for reflective BCs bottom of domain
   	//u(i, -2, k) = (k == YMOMENTUM ? -1.0 : 1.0) * u(i, 1, k);   
    //u(i, -1, k) = (k == YMOMENTUM ? -1.0 : 1.0) * u(i, 0, k);
				}
				//else if (i < 1 || i > (nozzlenumber * 2 - 2) * nozzleD * cellData  && i <= (nozzlenumber * 2 - 1) * nozzleD * cellData) {
			    else if (i >= (nozzlenumber - 1) * (wallD + nozzleD) *cellData && i <= ((nozzlenumber - 1) * (wallD + nozzleD) + wallD)*cellData - 1) {
							
					//u(i, -2, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 1, k);
					//u(i, -1, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 0, k);	
				u(i, -1, k) = u(i, 0, k);
				u(i, -2, k) = u(i, 1, k);
				u(i, -1, XMOMENTUM) = -u(i, 0, XMOMENTUM);
				u(i, -2, XMOMENTUM) = -u(i, 1, XMOMENTUM);
			    u(i, -1, YMOMENTUM) = -u(i, 0, YMOMENTUM);
				u(i, -2, YMOMENTUM) = -u(i, 1, YMOMENTUM);
				}
				}
      }
	  }	  
	  else {
		  //bot
		  real temp[NUMBER_VARIABLES];
			temp[DENSITY] = 1.0;
			temp[XVELOCITY] = 0;
			temp[YVELOCITY] = 0.0;
			temp[PRESSURE] = 1.0;
			temp[LAMBDA0] = 0.0;
			temp[LAMBDA1] = 0.0;
			primitiveToConservativeInPlace(temp);
			const real r = 0.05;
			for (int k = 0; k < NUMBER_VARIABLES; k++) {
				u(i, u.activeNy() + 1, k) = r * temp[k] + (1.0 - r) * u(i, u.activeNy() - 2, k);
				u(i, u.activeNy()    , k) = r * temp[k] + (1.0 - r) * u(i, u.activeNy() - 1, k);                   
			}
    //   u(i, u.activeNy() + 1, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, u.activeNy() - 2, k);
    //   u(i, u.activeNy()    , k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, u.activeNy() - 1, k);
    // for periodic BCs
    //    u(i, u.activeNy() + 1, k) = u(i, 1, k);
    //    u(i, u.activeNy()    , k) =  u(i, 0, k);
    // for reflective BCs top of domain
    //u(i, u.activeNy() + 1, k) = (k == YMOMENTUM ? -1.0 : 1.0) * u(i, u.activeNy() - 2, k);
    //u(i, u.activeNy()    , k) = (k == YMOMENTUM ? -1.0 : 1.0) * u(i, u.activeNy() - 1, k);       
      }
    }
  
  if (YDIR && u.exists(0, i)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (downstream) {
        // u(-2, i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(1, i, k);
        // u(-1, i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(0, i, k);
        //u(-2, i, k) = u(1, i, k);
        //u(-1, i, k) = u(0, i, k); 
		//for periodic left BCs
		u(-2, i, k) = u( u.activeNx() - 2, i, k);
        u(-1, i, k) = u( u.activeNx() - 1, i, k); 
      } else {
        // u(u.activeNx() + 1, i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(u.activeNx() - 2, i, k);
        // u(u.activeNx()    , i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(u.activeNx() - 1, i, k);
        //u(u.activeNx() + 1, i, k) = u(u.activeNx() - 2, i, k);
        //u(u.activeNx()    , i, k) = u(u.activeNx() - 1, i, k);
		//for periodic right BCs
		u(u.activeNx() + 1, i, k) = u(1, i, k);
        u(u.activeNx()    , i, k) = u(0, i, k);
      }
    }
  }
}

template<bool XDIR>
__global__ void setSpecialBoundaryConditionsKernel(Mesh<GPU>::type u) {
//  const bool YDIR = !XDIR;

 // const int ki = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();
 // const int kj = blockIdx.y * blockDim.y + threadIdx.y - u.ghostCells();
//  const int k = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();
/*
  if (XDIR && u.exists(k, 0) && k < u.i(XCORNER)) {
    const int j = u.j(YCORNER);
    for (int n = 0; n < 2; n++) {
      u(k, j + n + 1) = u(k, j - n);
      u(k, j + n + 1, YMOMENTUM) = -u(k, j + n + 1, YMOMENTUM);
    }
  }
  else if (!XDIR && u.exists(0, k) && k > u.j(YCORNER)) {
    const int i = u.i(XCORNER);
    for (int n = 0; n < 2; n++) {
      u(i - n - 1, k) = u(i + n, k);
      u(i - n - 1, k, XMOMENTUM) = -u(i - n - 1, k, XMOMENTUM);
    }
  }

  if (XDIR && u.exists(k, 0) && k > u.i(block_x1) && k < u.i(block_x2)) {

    const int j_y1 = u.j(block_y1);
    for (int n = 0; n < 2; n++) {
      u(k, j_y1 + n + 1) = u(k, j_y1 - n);
      u(k, j_y1 + n + 1, YMOMENTUM) = -u(k, j_y1 + n + 1, YMOMENTUM);
    }

    const int j_y2 = u.j(block_y2);
    for (int n = 0; n < 2; n++) {
      u(k, j_y2 - n - 1) = u(k, j_y2 + n);
      u(k, j_y2 - n - 1, YMOMENTUM) = -u(k, j_y2 - n - 1, YMOMENTUM);
    }

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(block_y1) && k < u.j(block_y2)) {

    const int i_x1 = u.i(block_x1);
    for (int n = 0; n < 2; n++) {
      u(i_x1 + n + 1, k) = u(i_x1 - n, k);
      u(i_x1 + n + 1, k, XMOMENTUM) = -u(i_x1 + n + 1, k, XMOMENTUM);
    }

    const int i_x2 = u.i(block_x2);
    for (int n = 0; n < 2; n++) {
      u(i_x2 - n - 1, k) = u(i_x2 + n, k);
      u(i_x2 - n - 1, k, XMOMENTUM) = -u(i_x2 - n - 1, k, XMOMENTUM);
    }

  }
*/ 
}

