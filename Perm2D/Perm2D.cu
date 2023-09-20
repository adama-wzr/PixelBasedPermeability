#include "Perm2D.cuh"

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
		std::cout << "Width (pixels) = " << width << " Height (pixels) = " << height << " Channel = " << channel << std::endl;
		std::cout << "Porosity = " << simInfo.porosity << std::endl;
	}

	// Define mesh related parameters

	simInfo.numCellsX = width*opts.MeshAmp;			// Simulation Grid width in number of cells
	simInfo.numCellsY = height*opts.MeshAmp;		// Simulation Grid height in number of cells
	simInfo.nElements = simInfo.numCellsY*simInfo.numCellsX;	// Number of elements (total)

	int *Grid = (int *)malloc(sizeof(int)*simInfo.numCellsX*simInfo.numCellsY);		// Array that will hold binary domain (solid vs fluid)

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

	return 0;
}