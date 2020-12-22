/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "ghost.hpp"
#include <sys/times.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>

bool halt = false;

void signalHandler(int signal = 0) {
  halt = true;
}

int main(int argc, char** argv) {
  // capture SIGINT
  signal(SIGINT, signalHandler);

  if (argc < 2) {
    std::cerr << "Invoke with " << argv[0] << " <configuration file>" << std::endl;
    exit(1);
  }

  // select a device
  int num_devices;
  cudaGetDeviceCount(&num_devices);
  std::cout << "#Found " << num_devices << " GPGPUs" << std::endl;
  cudaDeviceProp properties;
  int best_device = 0;
  if (num_devices > 1) {
    // if there's more than one, pick the one with the highest compute capability
    int best_computemode = 0, computemode;
    for (int device = 0; device < num_devices; device++) {
      cudaGetDeviceProperties(&properties, device);
      std::cout << "  #" << device << " " << properties.name << ": " << properties.multiProcessorCount << " processors, compute capability " << properties.major << "." << properties.minor << std::endl;
      computemode = properties.major << 4 + properties.minor;
      if (best_computemode < computemode) {
        best_computemode = computemode;
        best_device = device;
      }
    }
  }
        best_device = atoi(argv[2]);
  cudaGetDeviceProperties(&properties, best_device);
  std::cout << "#  using #" << best_device << " (" << properties.name << ")" << std::endl;
  cudaSetDevice(best_device);

  // start a timer to get the total wall time at the end of the run
  struct tms startTimes, endTimes;
  timespec startClock, endClock;
  times(&startTimes);
  clock_gettime(CLOCK_REALTIME, &startClock);

  Solver solver(argv[1]);

  Solver::Status status = Solver::OUTPUT;

  Mesh<CPU>::type uCPU(*solver.u, Mesh<CPU>::type::Allocate);

#ifdef GLOUTPUT
  OpenGLOutputter outputter(argc, argv, *solver.u);
  boost::thread outputterThread(boost::ref(outputter));
#endif

  // open the data file and output a header line
  std::stringstream filename;
  filename << solver.outputDirectory << "data";

  std::ofstream dataFile;
  dataFile.open(filename.str().c_str());

#ifdef GHOST
  dataFile << "#" << solver.u->time() << "\t\t";
  for (int i = 0; i < solver.geometries.size(); i++) {
    if (solver.geometries[i].rotating) {
      dataFile << "torque on level set " << i << "\t\t";
    }
  }
#endif
  for (int i = 0; i < solver.outputRadii.size(); i++) {
    dataFile << "P(r=" << solver.outputRadii[i] << ")\t\t";
    dataFile << "flux(r=" << solver.outputRadii[i] << ")\t\t";
  }
  dataFile << std::endl;

  do {
#ifdef GLOUTPUT
    if (solver.getStepNumber() % 1 == 0) {
      outputter.dispatchDraw(*solver.u);
      //outputter.paused = true;
      while (outputter.gridToRender != NULL); // SPIN
    }
#endif

    /*if (solver.getStepNumber() % 10 == 0) {
      dataFile << std::scientific << std::setprecision(10) << solver.u->time() << " ";
#ifdef GHOST
      for (int i = 0; i < solver.geometries.size(); i++) {
        if (solver.geometries[i].rotating) {
          dataFile  << solver.getTorque(solver.geometries[i]) << " ";
        }
      }
      for (int i = 0; i < solver.outputRadii.size(); i++) {
        std::pair<real, real> integrals = solver.getPressureIntegral(solver.outputRadii[i]);
        dataFile << integrals.first << " " << integrals.second << " ";
      }
#endif
      dataFile << std::endl;
    }*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Advance the frame
uCPU = *solver.u;

//int I_Check = (int) (X2+0.6*AMP) / uCPU.dy();

//Copy and paste///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 if (uCPU(1, 2)[PMAX] > (1.0+1.0e-4) && Check > 0.0){
  
  std::cout << "Copying data!";
  Check = -1.0;

// Read in density /////////////////////////////////////////////////
  std::ifstream inFile;
  inFile.open("Density_input000001.vtk");  

  for (int n_skip = 0; n_skip < Skip_lines; ++n_skip) {//ignore the subtitles
    inFile.ignore(256,'\n');
  }


    for (int nj = 0; nj < yinput * cellData; ++nj) {
//    for (int nj = cell*YCORNER; nj > 0; --nj) {
    for (int ni = 0 ; ni < xinput * cellData; ++ni) {

      inFile >> uCPU(ni, nj)[DENSITY];

      
    }
  }
 
  inFile.close();
 
// Read in x-momentum //////////////////////////////////////////////
  std::ifstream inFile1;
  inFile1.open("Xmom_input000001.vtk");  

  for (int n_skip = 0; n_skip < Skip_lines; ++n_skip) {
    inFile1.ignore(256,'\n');
  }
  
for (int nj = 0; nj < yinput * cellData; ++nj) {
	//for (int nj = 0; nj < uCPU.activeNy(); ++nj) {
//    for (int nj = cell*YCORNER1; nj < cell*YCORNER; ++nj) {
//    for (int nj = cell*YCORNER; nj > 0; --nj) {
    for (int ni = 0 ; ni < xinput * cellData; ++ni) {
      
      inFile1 >> uCPU(ni, nj)[XMOMENTUM];
      
    }

  }
 
  inFile1.close();

// Read in y-momentum //////////////////////////////////////////////
  std::ifstream inFile2;
  inFile2.open("Ymom_input000001.vtk");  

  for (int n_skip = 0; n_skip < Skip_lines; ++n_skip) {
    inFile2.ignore(256,'\n');
  }
  
for (int nj = 0; nj < yinput * cellData; ++nj) {
//    for (int nj = cell*YCORNER1; nj < cell*YCORNER; ++nj) {  
//    for (int nj = cell*YCORNER; nj > 0; --nj) {
    for (int ni = 0 ; ni < xinput * cellData; ++ni) {
      
      inFile2 >> uCPU(ni, nj)[YMOMENTUM];
     
    }
  }
  
  inFile2.close();

// Read in lambda_0 //////////////////////////////////////////////
  std::ifstream inFile3;
  inFile3.open("La0_input000001.vtk");  

  for (int n_skip = 0; n_skip < Skip_lines; ++n_skip) {
    inFile3.ignore(256,'\n');
  }
  
for (int nj = 0; nj < yinput * cellData; ++nj) {
//    for (int nj = cell*YCORNER1; nj < cell*YCORNER; ++nj) {  
//    for (int nj = cell*YCORNER; nj > 0; --nj) {
    for (int ni = 0 ; ni < xinput * cellData; ++ni) {
      
      inFile3 >> uCPU(ni, nj)[LAMBDA0];
      
    }
  }
  
  inFile3.close();

// Read in lambda_1 //////////////////////////////////////////////
  std::ifstream inFile4;
  inFile4.open("La1_input000001.vtk");  

  for (int n_skip = 0; n_skip < Skip_lines; ++n_skip) {
    inFile4.ignore(256,'\n');
  }
  
for (int nj = 0; nj < yinput * cellData; ++nj) {
//    for (int nj = cell*YCORNER1; nj < cell*YCORNER; ++nj) {  
//    for (int nj = cell*YCORNER; nj > 0; --nj) {
    for (int ni = 0 ; ni < xinput * cellData; ++ni) {
      
      inFile4 >> uCPU(ni, nj)[LAMBDA1];
      
    }
  }
  
  inFile4.close();

// Read in energy //////////////////////////////////////////////
  std::ifstream inFile5;
  inFile5.open("Energy_input000001.vtk");  

  for (int n_skip = 0; n_skip < Skip_lines; ++n_skip) {
    inFile5.ignore(256,'\n');
  }
  
for (int nj = 0; nj < yinput * cellData; ++nj) {
//    for (int nj = cell*YCORNER1; nj < cell*YCORNER; ++nj) {  
//    for (int nj = cell*YCORNER; nj > 0; --nj) {
    for (int ni = 0 ; ni < xinput * cellData; ++ni) {
      
      inFile5 >> uCPU(ni, nj)[ENERGY];
      
    }
  }
  
  inFile5.close();

// Read in Pmax//////////////////////////////////////////////
  std::ifstream inFile6;
  inFile6.open("Pmax_input000001.vtk");  

  for (int n_skip = 0; n_skip < Skip_lines; ++n_skip) {
    inFile6.ignore(256,'\n');
  }

for (int nj = 0; nj < yinput * cellData; ++nj) {
//    for (int nj = cell*YCORNER1; nj < cell*YCORNER; ++nj) {  
//    for (int nj = cell*YCORNER; nj > 0; --nj) {
    for (int ni = 0 ; ni < xinput * cellData; ++ni) {
      
      inFile6 >> uCPU(ni, nj)[PMAX];
     
    }
  }
  
  inFile6.close();
} 
 /*
 if (uCPU(uCPU.activeNx()-400, 200)[PMAX] > (1.0+1.0e-4)){  
  
  for (int ni = 0; ni < uCPU.activeNx()/3; ++ni) {
  for (int nj = 0; nj < uCPU.activeNy(); ++nj) {  
 
    uCPU(ni, nj)[DENSITY] = uCPU(ni+2*uCPU.activeNx()/3, nj)[DENSITY];
    uCPU(ni, nj)[XMOMENTUM]  = uCPU(ni+2*uCPU.activeNx()/3, nj)[XMOMENTUM];
    uCPU(ni, nj)[YMOMENTUM] = uCPU(ni+2*uCPU.activeNx()/3, nj)[YMOMENTUM];
    uCPU(ni, nj)[LAMBDA0] = uCPU(ni+2*uCPU.activeNx()/3, nj)[LAMBDA0];
    uCPU(ni, nj)[LAMBDA1] = uCPU(ni+2*uCPU.activeNx()/3, nj)[LAMBDA1];
    uCPU(ni, nj)[ENERGY]  = uCPU(ni+2*uCPU.activeNx()/3, nj)[ENERGY];
    uCPU(ni, nj)[PMAX] = uCPU(ni+2*uCPU.activeNx()/3, nj)[PMAX];

 }
 }
//Initialize new domain///////////////////////////////////////////////////////////////////////////////////////////// 

  for (int ni = uCPU.activeNx()/3; ni < uCPU.activeNx(); ++ni){
  for (int nj = 0; nj < uCPU.activeNy(); ++nj){

     uCPU(ni, nj)[DENSITY]    = 1.0;
     uCPU(ni, nj)[XMOMENTUM]  = 0.0;
     uCPU(ni, nj)[YMOMENTUM]  = 0.0;
     uCPU(ni, nj)[LAMBDA0]    = 0.0;
     uCPU(ni, nj)[LAMBDA1]    = 0.0;
     uCPU(ni, nj)[ENERGY]     = 1.0 / (gamma() - 1.0);
     uCPU(ni, nj)[PMAX]       = 1.0;                            
  }
  } 

//Add inert layer///////////////////////////////////////////////////////////////////////////////////////////////// 
  if (Counter > frame_interface ) {

	//	density_new = 1.0 + 0.1 * sin (2.0*PI*x/500);
     
		for (int ni = uCPU.activeNx()/3; ni < uCPU.activeNx(); ++ni){
		for (int nj = 0; nj < uCPU.activeNy(); ++nj){

		//uCPU(ni, nj)[LAMBDA0]    = 1.0;
		//uCPU(ni, nj)[LAMBDA1]    = 1.0;
		uCPU(ni, nj)[DENSITY]    = 1.0 + 0.25 * sin (2.0*(ni-uCPU.activeNx()/3)*PI/(cellData*wavelength));
		//uCPU(ni, nj)[DENSITY]    = 0.75;
                         
		}
		}

  }
                      

 Counter = Counter + 1.0;
} 
*/
*solver.u = uCPU;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (status == Solver::OUTPUT) {
      uCPU = *solver.u;

      if (true) {
/*                 
        std::stringstream filename;
        filename << solver.outputDirectory << "PRESSURE" << std::setw(6) << std::setfill('0') << solver.getOutputNumber() << ".vtk";

        std::ofstream outFile;
        outFile.open(filename.str().c_str());

        outFile.precision(8);
        outFile << "DIMENSIONS " << uCPU.activeNx() + 1 << " " << uCPU.activeNy() + 1<< " " << 1 << std::endl;
        outFile << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
        outFile << "SPACING " << uCPU.dx() << " " << uCPU.dy() << " 1" << std::endl;
        outFile << "CELL_DATA " << uCPU.activeNx() * uCPU.activeNy() << std::endl;

        outFile << "SCALARS pressure float" << std::endl;
        outFile << "LOOKUP_TABLE default" << std::endl;
      //  for (int j = 0; j < uCPU.activeNy(); ++j) { 
          int j = 1;
          for (int i = 0; i < uCPU.activeNx(); ++i) {
            real p[NUMBER_VARIABLES];
            conservativeToPrimitive(uCPU(i, j), p);
            outFile << std::fixed << p[PRESSURE] << std::endl;
          }
        // }
        outFile.close(); 
 */    
      
        std::stringstream filename;
        filename << solver.outputDirectory << "output" << std::setw(6) << std::setfill('0') << solver.getOutputNumber() << ".vtk";

        std::ofstream outFile;
        outFile.open(filename.str().c_str());

        outFile.precision(8);
 /*       outFile << "# vtk DataFile Version 3.0" << std::endl;
        outFile << "x y p Energy rho T" << std::endl;
        outFile << "ASCII" << std::endl;
        outFile << "DATASET STRUCTURED_POINTS" << std::endl;
        outFile << "DIMENSIONS " << uCPU.activeNx() + 1 << " " << uCPU.activeNy() + 1<< " " << 1 << std::endl;
        outFile << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
        outFile << "SPACING " << uCPU.dx() << " " << uCPU.dy() << " 1" << std::endl;
        outFile << "CELL_DATA " << uCPU.activeNx() * uCPU.activeNy() << std::endl;
*/
       /* outFile << "SCALARS density float" << std::endl;
        outFile << "LOOKUP_TABLE default" << std::endl;

        for (int j = 0; j < uCPU.activeNy(); ++j) {
          for (int i = 0; i < uCPU.activeNx(); ++i) {
            if (((uCPU.x(i) >= XCORNER2) || (uCPU.x(i) <= XCORNER1)) && (uCPU.y(j) >= YCORNER))
              outFile << 0 << std::endl;
            else
            outFile << std::fixed << uCPU(i, j)[DENSITY] << std::endl;
          }
        }
*/
       /* for (int j = 0; j < uCPU.activeNy(); ++j) {
          for (int i = 0; i < uCPU.activeNx(); ++i) {
            real p[NUMBER_VARIABLES];
            conservativeToPrimitive(uCPU(i, j), p);
            if (((uCPU.x(i) >= XCORNER2) || (uCPU.x(i) <= XCORNER1)) && (uCPU.y(j) >= YCORNER))
              outFile << 0 << std::endl;
            else
            outFile << std::fixed << p[PRESSURE] << std::endl;
          }
        }*/
	/*	double maxp=1.0;
		double maxt=1.0;
		for (int j = 0; j < 500; ++j) {
		for(int i=0;i < uCPU.activeNx(); ++i) {
			real p[NUMBER_VARIABLES];
            conservativeToPrimitive(uCPU(i, j), p);
			if (p[PRESSURE] > maxp) {
				maxp=p[PRESSURE];
				maxt=p[PRESSURE]/p[DENSITY];
			}
		}
		}
			outFile << std::fixed << maxp << " " << maxt << std::endl;
	*/		
	    int j=750;
		for(int i=0;i < uCPU.activeNx(); ++i) {
			real p[NUMBER_VARIABLES];
            conservativeToPrimitive(uCPU(i, j), p);
			outFile << std::fixed << i <<" "<< j <<" "<< p[PRESSURE] << " " << p[DENSITY] << " " << p[PRESSURE]/p[DENSITY] << std::endl;
		}			
		
/*
        outFile << "SCALARS soundSpeed float" << std::endl;
        outFile << "LOOKUP_TABLE default" << std::endl;
        for (int j = 0; j < uCPU.activeNy(); ++j) {
          for (int i = 0; i < uCPU.activeNx(); ++i) {
            if (((uCPU.x(i) >= XCORNER2) || (uCPU.x(i) <= XCORNER1)) && (uCPU.y(j) >= YCORNER))
              outFile << 0 << std::endl;
            else
            outFile << std::fixed << soundSpeed(uCPU(i, j)) << std::endl;
          }
        }

#ifdef GHOST
        outFile << "SCALARS SDF float" << std::endl;
        outFile << "LOOKUP_TABLE default" << std::endl;
        for (int j = 0; j < uCPU.activeNy(); ++j) {
          for (int i = 0; i < uCPU.activeNx(); ++i) {
            if (((uCPU.x(i) >= XCORNER2) || (uCPU.x(i) <= XCORNER1)) && (uCPU.y(j) >= YCORNER))
              outFile << 0 << std::endl;
            else
            outFile << std::fixed << uCPU(i, j)[PHI] << std::endl;
          }
        }
#endif

#ifdef REACTIVE
        outFile << "SCALARS fraction float" << std::endl;
        outFile << "LOOKUP_TABLE default" << std::endl;
        for (int j = 0; j < uCPU.activeNy(); ++j) {
          for (int i = 0; i < uCPU.activeNx(); ++i) {
            if (((uCPU.x(i) >= XCORNER2) || (uCPU.x(i) <= XCORNER1)) && (uCPU.y(j) >= YCORNER))
              outFile << 0 << std::endl;
            else
            outFile << std::fixed << uCPU(i, j)[Y] / uCPU(i, j)[DENSITY] << std::endl;
          }
        }

        outFile << "SCALARS shock float" << std::endl;
        outFile << "LOOKUP_TABLE default" << std::endl;
        for (int j = 0; j < uCPU.activeNy(); ++j) {
          for (int i = 0; i < uCPU.activeNx(); ++i) {
            if (((uCPU.x(i) >= XCORNER2) || (uCPU.x(i) <= XCORNER1)) && (uCPU.y(j) >= YCORNER))
              outFile << 0 << std::endl;
            else
            outFile << std::fixed << uCPU(i, j)[ISSHOCK] << std::endl;
          }
        }
#endif

        outFile << "VECTORS velocity float" << std::endl;
        for (int j = 0; j < uCPU.activeNy(); ++j) {
          for (int i = 0; i < uCPU.activeNx(); ++i) {
            real p[NUMBER_VARIABLES];
            conservativeToPrimitive(uCPU(i, j), p);
            if (((uCPU.x(i) >= XCORNER2) || (uCPU.x(i) <= XCORNER1)) && (uCPU.y(j) >= YCORNER))
              outFile << "0 0 0" << std::endl;
            else
            outFile << std::fixed << p[XVELOCITY] << " " << p[YVELOCITY] << " " << 0.0 << std::endl;
          }
        }*/
        outFile.close();
        
        std::stringstream filename1;
        filename1 << solver.outputDirectory << "MFR" << std::setw(6) << std::setfill('0') << solver.getOutputNumber() << ".vtk";

        std::ofstream outFile1;
        outFile1.open(filename1.str().c_str());

        outFile1.precision(8);
        //outFile1 << "# vtk DataFile Version 3.0" << std::endl;
        //outFile1 << "vtk output Mass Flow Rate" << std::endl;
        //outFile1 << "ASCII" << std::endl;
        //outFile1 << "DATASET STRUCTURED_POINTS" << std::endl;
        //outFile1 << "DIMENSIONS " << uCPU.activeNx() + 1 << " " << uCPU.activeNy() + 1<< " " << 1 << std::endl;
        //outFile1 << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
        //outFile1 << "SPACING " << uCPU.dx() << " " << uCPU.dy() << " 1" << std::endl;
        //outFile1 << "CELL_DATA " << uCPU.activeNx() * uCPU.activeNy() << std::endl;
		double MFR=0.0000000;		
		for (int nozzlenumber = 1; nozzlenumber <= slotNumber ; nozzlenumber++) {
				int j=0;
				
				for (int i = (((nozzlenumber - 1) * (wallD + nozzleD) + wallD)*cellData); i <= (nozzlenumber * (wallD + nozzleD) * cellData - 1); ++i) { 
				real p[NUMBER_VARIABLES];
                conservativeToPrimitive(uCPU(i, j), p);
				
				MFR = MFR + p[DENSITY] * sqrt(pow(p[XVELOCITY],2)+pow(p[YVELOCITY],2)) / (cellData);
				
				}
		}		
		outFile1 << std::fixed << MFR << std::endl;

		outFile1.close();
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		
		std::stringstream filename2;
        filename2 << solver.outputDirectory << "rhoV" << std::setw(6) << std::setfill('0') << solver.getOutputNumber() << ".vtk";

        std::ofstream outFile2;
        outFile2.open(filename2.str().c_str());

        outFile2.precision(8);
        //outFile1 << "# vtk DataFile Version 3.0" << std::endl;
        //outFile1 << "vtk output Mass Flow Rate" << std::endl;
        //outFile1 << "ASCII" << std::endl;
        //outFile1 << "DATASET STRUCTURED_POINTS" << std::endl;
        //outFile1 << "DIMENSIONS " << uCPU.activeNx() + 1 << " " << uCPU.activeNy() + 1<< " " << 1 << std::endl;
        //outFile1 << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
        //outFile1 << "SPACING " << uCPU.dx() << " " << uCPU.dy() << " 1" << std::endl;
        //outFile1 << "CELL_DATA " << uCPU.activeNx() * uCPU.activeNy() << std::endl;
		
		for (int nozzlenumber = 1; nozzlenumber <= slotNumber ; nozzlenumber++) {
				int j=0;
				
				for (int i = (((nozzlenumber - 1) * (wallD + nozzleD) + wallD)*cellData); i <= (nozzlenumber * (wallD + nozzleD) * cellData - 1); ++i) { 
				real p[NUMBER_VARIABLES];
                conservativeToPrimitive(uCPU(i, j), p);
				
				outFile2 << std::fixed << i <<" "<< j << " " << p[DENSITY] << " " << (sqrt(pow(p[XVELOCITY],2)+pow(p[YVELOCITY],2))) << std::endl;
				
				}
		}

		outFile2.close();
//////////////////////////////////////////////////////////////////////////////////////

/*	
       std::stringstream filename4;
        filename4 << solver.outputDirectory << "1D_Average_Profile" << std::setw(6) << std::setfill('0') << solver.getOutputNumber() << ".vtk";

        std::ofstream outFile4;
        outFile4.open(filename4.str().c_str());

        outFile4.precision(8);
       outFile4 << "DIMENSIONS " << uCPU.activeNx() + 1 << " " << uCPU.activeNy() + 1<< " " << 1 << std::endl;
        outFile4 << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
        outFile4 << "SPACING " << uCPU.dx() << " " << uCPU.dy() << " 1" << std::endl;
        outFile4 << "CELL_DATA " << uCPU.activeNx() * uCPU.activeNy() << std::endl;

        outFile4 << "1.P_avg, 2.rho_avg, 3.u_avg, 4.Z_avg" << std::endl;
   

        
        for (int i = 0; i < uCPU.activeNx(); ++i) {
          real RHO_SUM = 0.0, RHO_AVG = 0.0, P_SUM = 0.0, P_AVG = 0.0, U_RHO_SUM = 0.0; 
          real U_RHO_AVG = 0.0, Z_RHO_SUM = 0.0, Z_RHO_AVG = 0.0;
          for (int j = 0; j < uCPU.activeNy(); ++j) {
            real p[NUMBER_VARIABLES];
            conservativeToPrimitive(uCPU(i, j), p);
            RHO_SUM = RHO_SUM+uCPU(i, j)[DENSITY]*uCPU.dy();
            P_SUM = P_SUM+p[PRESSURE]*uCPU.dy();
            U_RHO_SUM = U_RHO_SUM+p[XVELOCITY]*uCPU(i, j)[DENSITY]*uCPU.dy();
            Z_RHO_SUM = Z_RHO_SUM+uCPU(i, j)[LAMBDA1]*uCPU.dy();            
          }
          P_AVG = P_SUM/(uCPU.activeNy()*uCPU.dy());
          RHO_AVG = RHO_SUM/(uCPU.activeNy()*uCPU.dy());
          U_RHO_AVG = U_RHO_SUM/(uCPU.activeNy()*uCPU.dy());
          Z_RHO_AVG = Z_RHO_SUM/(uCPU.activeNy()*uCPU.dy());
          outFile4 << std::fixed << P_AVG << "  " << RHO_AVG << "  " << U_RHO_AVG/RHO_AVG << "  " << Z_RHO_AVG/RHO_AVG << std::endl;
        }
     outFile4.close(); 
*/
//////////////////////////////////////////////////////////////////////////////////////
      }

      if (true) {
	  //if (Check < 0.0){

       if (Plot_counter > Plot_step) {
        Plot_counter = 1.0; 
        std::vector<ImageOutputs> outputs;

        outputs.push_back((ImageOutputs){"pressure", PRESSUREFOO, HUE, 1.0, 20.0});
        //outputs.push_back((ImageOutputs){"sootfoil", PMAXPLOT, CUBEHELIX, 10.0, 100.0});
       // outputs.push_back((ImageOutputs){"velocity_ y", XVELOCITYPLOT, HUE, -6.0, 6.0});
       //outputs.push_back((ImageOutputs){"pressure2", PRESSUREFOO, HUE, 1.0, 50.0});
        outputs.push_back((ImageOutputs){"density", DENSITY, HUE, 1.0, 10.0});
       // outputs.push_back((ImageOutputs){"LAMBDA0", LAMBDA0PLOT, HUE, 0.0, 1.0});
       outputs.push_back((ImageOutputs){"LAMBDA1", LAMBDA1PLOT, HUE, 0.0, 1.0});
       outputs.push_back((ImageOutputs){"schlieren", SCHLIEREN, GREYSCALE, 0.0, 1.0});
       outputs.push_back((ImageOutputs){"temperature", SOUNDSPEED, HUE, 4.0, 13.0});

        for (std::vector<ImageOutputs>::iterator iter = outputs.begin(); iter != outputs.end(); ++iter) {
          std::stringstream filename;
          filename << solver.outputDirectory << (*iter).prefix << std::setw(6) << std::setfill('0') << solver.getOutputNumber() << ".png";
          saveFrame(*solver.u, (*iter).plotVariable, (*iter).colourMode, filename.str().c_str(), (*iter).min, (*iter).max);
        }
      } else {
        Plot_counter = Plot_counter + 1.0;
      }

    // }
	 }
    }
  } while ((status = solver.step()) != Solver::FINISHED && !halt);

  times(&endTimes);
  clock_gettime(CLOCK_REALTIME, &endClock);
  const double wallTime = (endClock.tv_sec - startClock.tv_sec) + (endClock.tv_nsec - startClock.tv_nsec) * 1e-9;

  std::cout << "CPU time, wall= " << std::setprecision(2) << std::fixed << wallTime << "s, user=" << (endTimes.tms_utime - endTimes.tms_utime) << "s, sys=" << (endTimes.tms_stime - endTimes.tms_stime) << "s.  Time for {fluxes=" << solver.getTimeFluxes() * 1e-3 << "s, sources=" << solver.getTimeSourcing() * 1e-3 << "s, reduction=" << solver.getTimeReducing() * 1e-3 << "s, adding=" << solver.getTimeAdding() * 1e-3 << "s}" << std::endl;
  return 0;
}


