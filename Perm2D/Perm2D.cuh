#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <stdbool.h>
#include <fstream>
#include <cfloat>
#include <set>
#include <string>


// Structure for saving user entered options
typedef struct
{
	float density;
	float viscosity;
	float PL;
	float PR;
	int MeshAmp;
	char* inputFilename;
	char* outputFilename;
	bool printMaps;
	bool verbose;
	long int MaxIterSolver;
	float ConvergenceSolver;
	long int MaxIterGlobal;
	float ConvergenceRMS;
	float DomainHeight;
	float DomainWidth;
}options;

// Struct to hold constants, variables, and results intrinsic to the simulation
typedef struct{
	int numCellsX;
	int numCellsY;
	int nElements;
	float porosity;
	float gpuTime;
	float Perm;
	float Convergence;
}simulationInfo;

int printOptions(options* opts){
	/*
		printOptions Function:
		Inputs:
			- Pointer to struct options: struct with the user entered options.
		Outputs: None

		Function prints all the user entered parameters.
	*/
	printf("--------------------------------------\n\n");
	printf("User Selected Options:\n\n");
	printf("--------------------------------------\n");
	printf("Width = %.2f m\n", opts->DomainWidth);
	printf("Height = %.2f m\n", opts->DomainHeight);
	printf("Pressure Left = %.2f Pa\n", opts->PL);
	printf("Pressure Right = %.2f Pa\n", opts->PR);
	printf("Density = %.2f kg/m^3\n", opts->density);
	printf("Kinematic Viscosity = %.2f m^2/s\n", opts->viscosity);
	printf("Mesh Refinement = %d\n", opts->MeshAmp);
	printf("Maximum Iterations Solver = %ld\n", opts->MaxIterSolver);
	printf("Solver Convergence = %.10f\n", opts->ConvergenceSolver);
	printf("Maximum Iterations Global = %ld\n", opts->MaxIterGlobal);
	printf("RMS Convergence = %.10f\n", opts->ConvergenceRMS);
	printf("Name of input image: %s\n", opts->inputFilename);
	printf("Name of output file: %s\n", opts->outputFilename);

	return 0;
}


int readInputFile(char* FileName, options* opts){

	/*
		readInputFile Function:
		Inputs:
			- FileName: pointer to where the input file name is stored.
			- Pointer to struct options: struct with the user entered options.
		Outputs: None

		Function reads the input file and stores the options in the opts struct.
	*/

	std::string myText;

	char tempC[1000];
	float tempD;
	char tempFilenames[1000];
	std::ifstream InputFile(FileName);

	// initialize the pointers so they are not random

	opts->inputFilename=(char*)malloc(100*sizeof(char));
	opts->outputFilename=(char*)malloc(100*sizeof(char));
	while(std::getline(InputFile, myText)){

	 	sscanf(myText.c_str(), "%s %f", tempC, &tempD);
	 	if (strcmp(tempC, "Dens:") == 0){
	 		opts->density = tempD;
	 	}else if(strcmp(tempC, "Visc:") == 0){
	 		opts->viscosity = tempD;

	 	}else if(strcmp(tempC, "MeshAmp:") == 0){
	 		opts->MeshAmp = (int)tempD;

	 	}else if(strcmp(tempC, "PL:") == 0){
	 		opts->PL = tempD;

	 	} else if(strcmp(tempC, "PR:") == 0){
	 		opts->PR = tempD;

	 	}else if(strcmp(tempC, "InputName:") == 0){
	 		sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
	 		strcpy(opts->inputFilename, tempFilenames);

	 	}else if(strcmp(tempC, "OutputName:") == 0){
	 		sscanf(myText.c_str(), "%s %s", tempC, tempFilenames);
	 		strcpy(opts->outputFilename, tempFilenames);

	 	}else if(strcmp(tempC, "printMaps:") == 0){
	 		opts->printMaps = (bool)tempD;

	 	}else if(strcmp(tempC, "MaxIterGlobal:") == 0){
	 		opts->MaxIterGlobal = (long int)tempD;

	 	}else if(strcmp(tempC, "ResidualConv:") == 0){
	 		opts->ConvergenceRMS = tempD;

	 	}else if(strcmp(tempC, "MaxIterSolver:") == 0){
	 		opts->MaxIterSolver = (long int)tempD;

	 	}else if(strcmp(tempC, "SolverConv:") == 0){
	 		opts->ConvergenceSolver = tempD;

	 	}else if(strcmp(tempC, "Verbose:") == 0){
	 		opts->verbose = (bool)tempD;
	 	}else if(strcmp(tempC, "DomainWidth:") == 0){
	 		opts->DomainWidth = tempD;

	 	}else if(strcmp(tempC, "DomainHeight:") == 0){
	 		opts->DomainHeight = tempD;

	 	}
	}
	
	InputFile.close();



	if(opts->verbose == 1){
		printOptions(opts);
	}
	return 0;
}

int readImage(unsigned char** imageAddress, int* Width, int* Height, int* NumOfChannels, char* ImageName){
	/*
		readImage Function:
		Inputs:
			- imageAddress: unsigned char reference to the pointer in which the image will be read to.
			- Width: pointer to variable to store image width
			- Height: pointer to variable to store image height
			- NumofChannels: pointer to variable to store number of channels in the image.
					-> NumofChannels has to be 1, otherwise code is terminated. Please enter grayscale
						images with NumofChannels = 1.
		Outputs: None

		Function reads the image into the pointer to the array to store it.
	*/

	*imageAddress = stbi_load(ImageName, Width, Height, NumOfChannels, 1);

	return 0;
}


float calcPorosity(unsigned char* imageAddress, int Width, int Height){
	/*
		calcPorosity
		Inputs:
			- imageAddress: pointer to the read image.
			- Width: original width from std_image
			- Height: original height from std_image

		Output:
			- porosity: float containing porosity.

		Function calculates porosity by counting pixels.
	*/

	float totalCells = (float)Height*Width;
	float porosity = 0;
	for(int i = 0; i<Height; i++){
		for(int j = 0; j<Width; j++){
			if(imageAddress[i*Width + j] < 150){
				porosity += 1.0/totalCells;
			}
		}
	}

	return porosity;
}