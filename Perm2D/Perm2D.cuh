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
}options;


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
	printf("Pressure Left = %.2f\n", opts->PL);
	printf("Pressure Right = %.2f\n", opts->PR);
	printf("Density = %.2f\n", opts->density);
	printf("Viscosity = %.2f\n", opts->viscosity);
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
	 	}
	}
	
	InputFile.close();



	if(opts->verbose == 1){
		printOptions(opts);
	}
	return 0;
}

