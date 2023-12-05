#include "Perm2D.cuh"

int main(void){

	// More efficient printing with parallel computing/Linux
	fflush(stdout);

	// Parse user entered options

	options opts;	// struct to hold options

	char inputFilename[100];

	sprintf(inputFilename, "input.txt");

	readInputFile(inputFilename, &opts);

	// Set OpenMP CPU environment

	int numThreads = opts.nCores;
	int nGPUs;

	omp_set_num_threads(numThreads);
	cudaGetDeviceCount(&nGPUs);

	printf("Do I even get here?CPU = %d, GPU = %d\n", numThreads, nGPUs);

	// Start datastructures

	options myOpts;

	simulationInfo simInfo;

	domainInfo info;

	#pragma omp parallel for schedule(runtime) private(myOpts, simInfo, info)
	for(int myImg = 0; myImg<opts.nImg; myImg++){

		// read options
		readInputFile(inputFilename, &myOpts);

		printf("Pressure Left = %f\n", myOpts.PL);
		
		// Get thread index

		int threadIdx = omp_get_thread_num();

		printf("Thread idx = %d, img num = %d\n", threadIdx, myImg);

		cudaSetDevice(threadIdx);

		// Read 2D Input Image

		int width, height, channel;
		unsigned char* targetImage;
		char filename[100];

		sprintf(filename, "%05d.jpg", myImg);

		readImage(&targetImage, &width, &height, &channel, filename);

		myOpts.inputFilename = filename;

		if (channel != 1){
			printf("Error: please enter a grascale image with 1 channel.\n Current number of channels = %d\n", channel);
		}

		simInfo.porosity = calcPorosity(targetImage, width, height);

		if(myOpts.verbose == 1){
			std::cout << "Image Parameters:" << std::endl;
			std::cout <<  "\n--------------------------------------" << std::endl;
			std::cout << "Width (pixels) = " << width << " Height (pixels) = " << height << " Channel = " << channel << std::endl;
			std::cout << "Porosity = " << simInfo.porosity << std::endl;
		}

		// Define mesh related parameters

		simInfo.numCellsX = width*myOpts.MeshAmp;			// Simulation Grid width in number of cells
		simInfo.numCellsY = height*myOpts.MeshAmp;		// Simulation Grid height in number of cells
		simInfo.nElements = simInfo.numCellsY*simInfo.numCellsX;	// Number of elements (total)
		simInfo.dx = myOpts.DomainWidth/simInfo.numCellsX;			// dx
		simInfo.dy = myOpts.DomainWidth/simInfo.numCellsY;			// dy

		unsigned int *Grid = (unsigned int*)malloc(sizeof(int)*simInfo.numCellsX*simInfo.numCellsY);		// Array that will hold binary domain (solid vs fluid)

		// Mesh Amplify and decode image into binary matrix

		for(int i = 0; i<simInfo.numCellsY; i++){
			for (int j = 0; j<simInfo.numCellsX; j++){
				int targetIndex_Row = i/myOpts.MeshAmp;
				int targetIndex_Col = j/myOpts.MeshAmp;
				if(targetImage[targetIndex_Row*width + targetIndex_Col] < 150){
					Grid[i*simInfo.numCellsX + j] = 0; 			// black => fluid => 0 => void
				} else{
					Grid[i*simInfo.numCellsX + j] = 1;			// white => solid => 1 => material
				}
			}
		}

		// Next step is running the pathfinding algorithm
		// If there is no path, don't simulate permeability

		bool pathFlag = false;

		
		info.xSize = simInfo.numCellsX;
		info.ySize = simInfo.numCellsY;
		info.verbose = 0;

		pathFlag = aStarMain(Grid, info);

		if(pathFlag == false){
			std::cout << "No valid path found, exiting now." << std::endl;
		}else if(myOpts.verbose == 1){
			std::cout << "Valid path found.\nProceeding to permeability CFD simulation." << std::endl;
		}

		// Flood Fill to eliminate non-participating media

		FloodFill(Grid, &simInfo);

		std::cout << "Flood Fill Successfull." << std::endl;

		// Define arrays essential for the solution

		float *Pressure = (float *)malloc(sizeof(float)*simInfo.nElements);					// store pressure
		float *U = (float *)malloc(sizeof(float)*(simInfo.numCellsX+1)*simInfo.numCellsY);	// store U velocity
		float *V = (float *)malloc(sizeof(float)*(simInfo.numCellsY+1)*simInfo.numCellsX);	// store V velocity

		float *uExp = (float *)malloc(sizeof(float)*(simInfo.numCellsX+1)*simInfo.numCellsY);	// store explicit U velocity
		float *vExp = (float *)malloc(sizeof(float)*(simInfo.numCellsY+1)*simInfo.numCellsX);	// store explicit V velocity

		float *uCoeff = (float *)malloc(sizeof(float)*(simInfo.numCellsX+1)*simInfo.numCellsY);	// store U velocity coefficients
		float *vCoeff = (float *)malloc(sizeof(float)*(simInfo.numCellsY+1)*simInfo.numCellsX);	// store V velocity coefficients

		std::cout << "Allocated arrays Successfull." << std::endl;

		// Initialize arrays

		for(int row = 0; row<simInfo.numCellsY; row++){
			for(int col = 0; col< simInfo.numCellsX+1; col++){
				int index = row*(simInfo.numCellsX + 1) + col;
				if(col < simInfo.numCellsX){
					// Pressure[row*(simInfo.numCellsX) + col] = (opts.PL + opts.PR)/2;
					Pressure[row*(simInfo.numCellsX) + col] =  (1.0 - (float)col/(simInfo.numCellsX))*(myOpts.PL - myOpts.PR) + myOpts.PR;
				}
				U[index] = 0.01;
				V[index] = 0.0;
				uExp[index] = 0.01;
				vExp[index] = 0.0;
				uCoeff[index] = 0.0;
				vCoeff[index] = 0.0;
			}	
		}

		// Now we use the SUV-CUT algorithm to solve velocity-pressure coupled

		float RMS = 1.0;
		long int iter = 0;

		std::cout << "Start loop Successfull." << std::endl;

		while(iter < myOpts.MaxIterGlobal && RMS > myOpts.ConvergenceRMS){

			/*
				SUV-CUT procedure:
				- Solve for explicit component of u and v velocities
				- Use explicit u and v to solve for pressure implicitly
				- Use pressure solutions to correct u and v explicitly

				- (optional) solve equations of state to update physical properties

				Repeat until converged.

			*/

			if(iter == 0){
				printf("Global Iter: %ld, Thread = %d\n\n", iter+1, threadIdx);
			}else if(iter % 10 == 0){
				printf("Global Iter: %ld, Thread = %d\n", iter+1, threadIdx);
				printf("Permeability: %f\n", simInfo.Perm);
				printf("Continuity RMS: %1.9f\n\n", RMS);
			}
			explicitMomentum(Grid, uExp, vExp, U, V, uCoeff, vCoeff, &myOpts, &simInfo);
			implicitPressure(Grid, uExp, vExp, uCoeff, vCoeff, Pressure, &myOpts, &simInfo);
			momentumCorrection(Grid, uExp, vExp, U, V, uCoeff, vCoeff, Pressure, &myOpts, &simInfo);

			RMS = ResidualContinuity(U, V, &myOpts, &simInfo);

			PermCalc(U, &myOpts, &simInfo);

			iter++;
		}

		if(myOpts.printMaps == 1){
			printPUVmaps(Pressure, U, V, &myOpts, &simInfo);
		}

		printBatchOut(&myOpts, &simInfo, myImg, iter, RMS);

		// Save results to output file
		
		// Housekeeping

		free(Pressure);
		free(U);
		free(V);
		free(uCoeff);
		free(vCoeff);
		free(uExp);
		free(vExp);
		free(Grid);
		free(targetImage);

		// reset simInfo and info structs

		memset(&simInfo, 0, sizeof(simInfo));
		memset(&info, 0, sizeof(info));

	}

	return 0;
}