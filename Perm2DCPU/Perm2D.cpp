#include "Perm2D.h"

int main(void){

	// More efficient printing with parallel computing/Linux
	fflush(stdout);

	// Parse user entered options

	options opts;	// struct to hold options
	simulationInfo simInfo;		// Struct to hold constants, variables, and results intrinsic to the simulation

	char inputFilename[100];

	sprintf(inputFilename, "input.txt");

	readInputFile(inputFilename, &opts);

	// Read 2D Input Image

	int width, height, channel;
	unsigned char* targetImage;

	readImage(&targetImage, &width, &height, &channel, opts.inputFilename);

	if (channel != 1){
		printf("Error: please enter a grascale image with 1 channel.\n Current number of channels = %d\n", channel);
		return 1;
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
	simInfo.dx = simInfo.numCellsX/opts.DomainWidth;			// dx
	simInfo.dy = simInfo.numCellsY/opts.DomainHeight;			// dy

	unsigned int *Grid = (unsigned int*)malloc(sizeof(int)*simInfo.numCellsX*simInfo.numCellsY);		// Array that will hold binary domain (solid vs fluid)

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

	// Next step is running the pathfinding algorithm
	// If there is no path, don't simulate permeability

	bool pathFlag = false;

	domainInfo info;
	info.xSize = simInfo.numCellsX;
	info.ySize = simInfo.numCellsY;
	info.verbose = 0;

	pathFlag = aStarMain(Grid, info);

	if(pathFlag == false){
		std::cout << "No valid path found, exiting now." << std::endl;
		return 1;
	}else if(opts.verbose == 1){
		std::cout << "Valid path found.\nProceeding to permeability CFD simulation." << std::endl;
	}

	// Define arrays essential for the solution

	float *Pressure = (float *)malloc(sizeof(float)*simInfo.nElements);					// store pressure
	float *U = (float *)malloc(sizeof(float)*(simInfo.numCellsX+1)*simInfo.numCellsY);	// store U velocity
	float *V = (float *)malloc(sizeof(float)*(simInfo.numCellsY+1)*simInfo.numCellsX);	// store V velocity

	float *uExp = (float *)malloc(sizeof(float)*(simInfo.numCellsX+1)*simInfo.numCellsY);	// store explicit U velocity
	float *vExp = (float *)malloc(sizeof(float)*(simInfo.numCellsY+1)*simInfo.numCellsX);	// store explicit V velocity

	float *uCoeff = (float *)malloc(sizeof(float)*(simInfo.numCellsX+1)*simInfo.numCellsY);	// store U velocity coefficients
	float *vCoeff = (float *)malloc(sizeof(float)*(simInfo.numCellsY+1)*simInfo.numCellsX);	// store V velocity coefficients

	// Initialize arrays

	memset(Pressure, (opts.PL + opts.PR)/2, sizeof(Pressure));		// Initialized to avg. between PL and PR

	memset(uExp, 0, sizeof(uExp));									// Initialized to 0 because we will solve for it first step
	memset(vExp, 0, sizeof(vExp));									// Initialized to 0 because we will solve for it first step

	memset(uCoeff, 0, sizeof(uCoeff));								// Initialized to 0 because we solve for it first step
	memset(vCoeff, 0, sizeof(vCoeff));								// Initialized to 0 because we solve for it first step

	memset(V, 0, sizeof(V));										// Initialize to 0 because it is the dominant flow
	memset(U, 0.001, sizeof(U));									// Initialized to 0.001 just to help convergence

	// Now we use the SUV-CUT algorithm to solve velocity-pressure coupled

	float RMS = 1;
	long int iter = 0;

	while(iter < opts.MaxIterGlobal && RMS > opts.ConvergenceRMS){

		/*
			SUV-CUT procedure:
			- Solve for explicit component of u and v velocities
			- Use explicit u and v to solve for pressure implicitly
			- Use pressure solutions to correct u and v explicitly

			Repeat until converged.

		*/

		explicitMomentum(Grid, uExp, vExp, U, V, uCoeff, vCoeff, &opts, &simInfo);

	}

	return 0;
}