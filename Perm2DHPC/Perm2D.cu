#include "Perm2D.cuh"

int main(void){

	// More efficient printing with parallel computing/Linux
	fflush(stdout);

	// Parse user entered options

	bool convFlag = true;

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

	// Let's pre allocate all arrays here globally and then distribute them

	int nRows = opts.MeshAmp*128;
	int nCols = opts.MeshAmp*128;
	int nElements = nRows*nCols;

	unsigned int *Global_Grid = (unsigned int*)malloc(sizeof(int)*nElements*nGPUs);

	float *Global_Pressure = (float *)malloc(sizeof(float)*nElements*nGPUs);					// store pressure
	float *Global_U = (float *)malloc(sizeof(float)*(nCols+1)*nRows*nGPUs);	// store U velocity
	float *Global_V = (float *)malloc(sizeof(float)*(nCols+1)*nRows*nGPUs);	// store V velocity

	float *Global_uExp = (float *)malloc(sizeof(float)*(nCols+1)*nRows*nGPUs);	// store explicit U velocity
	float *Global_vExp = (float *)malloc(sizeof(float)*(nCols+1)*nRows*nGPUs);	// store explicit V velocity

	float *Global_uCoeff = (float *)malloc(sizeof(float)*(nCols+1)*nRows*nGPUs);	// store U velocity coefficients
	float *Global_vCoeff = (float *)malloc(sizeof(float)*(nCols+1)*nRows*nGPUs);	// store V velocity coefficients

	#pragma omp parallel for schedule(auto)
	for(int myImg = 0; myImg<opts.nImg; myImg++){

		// Start datastructures

		simulationInfo simInfo;
		convInfo *Conv;

		Conv = (convInfo *)malloc(sizeof(convInfo)*opts.MaxIterGlobal);

		memset(Conv, 0, sizeof(convInfo)*opts.MaxIterGlobal);
		
		// Get thread index

		int threadIdx = omp_get_thread_num();

		printf("Thread idx = %d, Img num = %d\n", threadIdx, myImg);

		cudaSetDevice(threadIdx);

		// Read 2D Input Image

		int width, height, channel;
		unsigned char* targetImage;
		char filename[100];
		char convFile[100];

		sprintf(filename, "%05d.jpg", myImg);

		readImage(&targetImage, &width, &height, &channel, filename);

		if (channel != 1){
			printf("Error: please enter a grascale image with 1 channel.\n Current number of channels = %d\n", channel);
		}

		simInfo.porosity = calcPorosity(targetImage, width, height);

		if(opts.verbose == 1){
			std::cout << "Image Parameters:" << std::endl;
			std::cout <<  "\n--------------------------------------" << std::endl;
			std::cout << "Width (pixels) = " << width << " Height (pixels) = " << height << " Channel = " << channel << std::endl;
			std::cout << "Porosity = " << simInfo.porosity << std::endl;
		}

		// Define mesh related parameters

		simInfo.numCellsX = width*opts.MeshAmp;			// Simulation Grid width in number of cells
		simInfo.numCellsY = height*opts.MeshAmp;		// Simulation Grid height in number of cells
		simInfo.nElements = simInfo.numCellsY*simInfo.numCellsX;	// Number of elements (total)
		simInfo.dx = opts.DomainWidth/simInfo.numCellsX;			// dx
		simInfo.dy = opts.DomainWidth/simInfo.numCellsY;			// dy

		unsigned int *Grid = Global_Grid + threadIdx * nElements;		// Array that will hold binary domain (solid vs fluid)

		// Mesh Amplify and decode image into binary matrix

		for(int i = 0; i<simInfo.numCellsY; i++){
			for (int j = 0; j<simInfo.numCellsX; j++){
				int targetIndex_Row = i/opts.MeshAmp;
				int targetIndex_Col = j/opts.MeshAmp;
				if(targetImage[targetIndex_Row*width + targetIndex_Col] < 150){
					Grid[i*simInfo.numCellsX + j] = 0; 			// black => fluid => 0 => void
				} else{
					Grid[i*simInfo.numCellsX + j] = 1;			// white => solid => 1 => material
				}
			}
		}

		// Free target_image since it is not needed anymore
		free(targetImage);

		// Flood Fill to eliminate non-participating media

		FloodFill(Grid, &simInfo);

		std::cout << "Flood Fill Successfull. Thread = " << threadIdx << std::endl;

		// Define arrays essential for the solution

		float *Pressure = Global_Pressure + threadIdx * nElements;
		float *U = Global_U + threadIdx * (nCols + 1)*nRows;
		float *V = Global_V + threadIdx * (nRows + 1)*nCols;

		float *uExp = Global_uExp + threadIdx * (nCols + 1)*nRows;
		float *vExp = Global_vExp + threadIdx * (nRows + 1)*nCols;

		float *uCoeff = Global_uCoeff + threadIdx * (nCols + 1)*nRows;
		float *vCoeff = Global_vCoeff + threadIdx * (nRows + 1)*nCols;

		std::cout << "Allocated arrays Successfull. Thread = " << threadIdx << std::endl;

		// Initialize arrays

		for(int row = 0; row<simInfo.numCellsY; row++){
			for(int col = 0; col< simInfo.numCellsX+1; col++){
				int index = row*(simInfo.numCellsX + 1) + col;
				if(col < simInfo.numCellsX){
					Pressure[row*(simInfo.numCellsX) + col] =  (1.0 - (float)col/(simInfo.numCellsX))*(opts.PL - opts.PR) + opts.PR;
				}
				U[index] = 0.01;
				V[index] = 0.0;
				uExp[index] = 0.01;
				vExp[index] = 0.0;
				uCoeff[index] = 0.0;
				vCoeff[index] = 0.0;
			}	
		}

		// start file to get convergence data

		if(convFlag == true)
		{
			sprintf(convFile, "ConvData_%05d.csv", myImg);
		}

		// Now we use the SUV-CUT algorithm to solve velocity-pressure coupled

		float RMS = 1.0;
		long int iter = 0;

		float PermTHR = 0.001;
		float PermOld = 1;
		float PermChange = 1;

		std::cout << "Start loop Successfull." << std::endl;

		while(iter < opts.MaxIterGlobal && RMS > opts.ConvergenceRMS && PermChange > PermTHR){

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
			// std::cout << "Thread Num:" << omp_get_thread_num() << "Explicit Momentum Iter" << iter << std::endl;
			explicitMomentum(Grid, uExp, vExp, U, V, uCoeff, vCoeff, &opts, &simInfo);
			// std::cout << "Thread Num:" << omp_get_thread_num() << "implicitPressure Iter" << iter << std::endl;
			implicitPressure(Grid, uExp, vExp, uCoeff, vCoeff, Pressure, &opts, &simInfo);
			// std::cout << "Thread Num:" << omp_get_thread_num() << "Momentum Correction Iter" << iter << std::endl;
			momentumCorrection(Grid, uExp, vExp, U, V, uCoeff, vCoeff, Pressure, &opts, &simInfo);

			RMS = ResidualContinuity(U, V, &opts, &simInfo);

			PermCalc(U, &opts, &simInfo);

			// Update our convergence file

			Conv[iter].iter = iter;
			Conv[iter].Perm = simInfo.Perm;
			Conv[iter].Residual = RMS;
			Conv[iter].PermChange = PermChange;


			// Calculate Perm change to see if it flatlined over 100 iterations
			if(iter % 100 == 0){
				PermChange = fabs((simInfo.Perm - PermOld)/simInfo.Perm);
				PermOld = simInfo.Perm;
			}

			iter++;
		}

		// Print PUV Map

		if(opts.printMaps == 1){
			printPUVmaps(Pressure, U, V, &opts, &simInfo, myImg);
		}

		// Print batch output

		printBatchOut(&opts, &simInfo, myImg, iter, RMS);

		// Print convergence data if user wants it
		if(convFlag == true){
			FILE *CONV;
			CONV = fopen(convFile, "w+");
			fprintf(CONV, "iter,K,R,Kchange\n");
			for(int i=0; i<opts.MaxIterGlobal; i++){
				fprintf(CONV,"%ld,%f,%f,%f\n", Conv[i].iter, Conv[i].Perm, Conv[i].Residual, Conv[i].PermChange);
			}
			fclose(CONV);
		}

		free(Conv);

	}

	return 0;
}