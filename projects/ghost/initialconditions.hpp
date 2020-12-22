/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#pragma once

__global__ void setInitialConditions(Mesh<GPU>::type u) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x,
            j = blockIdx.y * blockDim.y + threadIdx.y;

  if (u.active(i, j)) {
    real lambda0 = 0, lambda1 = 0, density = 0, v_x = 0, v_y = 0, p = 0, FLOOR_real_x = 0.0, FLOOR_real_y = 0.0;
    const real x = u.x(i), y = u.y(j);
    
   // if (x < 150 &&  y < 192) {
    if (x < 150 && y < 192) {
  //     if (x > X1 && x < X2) {
       p = pShock;
       density = rhoShock;
	   
	//   v_x = sqrt(1.32)*Mcj;
//v_x=2*1.053805645845247;
	   //lambda0 = 1;
	   
    }  
		/* if(x < 40.0 && (x-40)*(x-40) + y*y < 40*40) //CIRCLE
	{
	p = 2*pShock;
	density = 2*rhoShock;
	v_x = 2*1.053805645845247;
//	lambda1=0;
	} */
	/*
	else if (x > 400 && x < 1200) {
	density = 1.0 + 0.25 * sin (2.0*(x-400)*PI/wavelength);
	//density = 0.75;
	p=1.0;
	}
	*/
	/*
    else if (x > 192 && x < 2.0 * sin (2.0*PI*y/20)+202 && y > 0 && y < 192) {
  //     if (x > X1 && x < X2) {
       p = 1.0;
       density = 1.0 + 0.25 * sin (2.0*x*PI/10);
	   
	   //v_x = Mcj;
	   //lambda0 = 1;
	   
    } 
	*/
	else {
		p = 1.0;
		density = 1.0;			
	}			

    u(i, j, DENSITY)   = density;
    u(i, j, XMOMENTUM) = v_x * density;
    u(i, j, YMOMENTUM) = v_y * density;
    u(i, j, LAMBDA0)   = lambda0 * density;
    u(i, j, LAMBDA1)   = lambda1 * density;
    u(i, j, ENERGY)    = p / (gamma() - 1.0) + 0.5 * density * (v_x * v_x + v_y * v_y);
    u(i, j, PMAX) = p;
  }
}

