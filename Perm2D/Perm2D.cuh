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
	float dx;
	float dy;
	float porosity;
	float gpuTime;
	float Perm;
	float Convergence;
}simulationInfo;


//All data structures below are for a* algorithm:

// generate structure to store global information about the domain

typedef struct{
	unsigned int xSize;
	unsigned int ySize;
	bool verbose;
}domainInfo;

// node struct will hold information on cell parent cells, f, g, and h.

typedef struct{
	int parentRow, parentCol;
	float f,g,h;
}node;

// define pair for coords

typedef std::pair<int, int> coordPair;

// define pair <float, pair<i,j>> for open list

typedef std::pair<float, std::pair<int,int> > OpenListInfo;

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
	printf("--------------------------------------\n\n");

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


int aStarMain(unsigned int* GRID, domainInfo info){
	/*
		aStarMain Function
		Inputs:
			- unsigned int GRID: grid, at each location either a 1 or a 0.
				1 means solid, 0 void. Those are the boundary conditions.
			- domainInfo info: data structure with info regarding height and width of the domain
		Output:
			- either a one or a zero. One means there is a path, zero means there isn't.
	*/

	// Initialize both lists, open and closed as arrays

	bool* closedList = (bool *)malloc(sizeof(bool)*info.xSize*info.ySize);

	memset(closedList, false, sizeof(closedList));

	// Declare 2D array of structure type "node"
	// Node contains information such as parent coordinates, g, h, and f

	node nodeInfo[info.ySize][info.ySize];

	// Initialize all paremeters

	for(int i = 0; i<info.ySize; i++){
		for(int j = 0; j<info.xSize; j++){
			nodeInfo[i][j].f = FLT_MAX;
			nodeInfo[i][j].g = FLT_MAX;
			nodeInfo[i][j].h = FLT_MAX;
			nodeInfo[i][j].parentCol = -1;
			nodeInfo[i][j].parentRow = -1;
		}
	}

	// Initialize parameters for all starting nodes

	for(int i = 0; i<info.ySize; i++){
		if(GRID[i*info.xSize + 0] == 0){
			nodeInfo[i][0].f = 0.0;
			nodeInfo[i][0].g = 0.0;
			nodeInfo[i][0].h = 0.0;
			nodeInfo[i][0].parentCol = 0;
			nodeInfo[i][0].parentRow = i;
		}
	}

	// Create open list

	std::set<OpenListInfo> openList;

	// Insert all starting nodes into the open list

	for(int i = 0; i<info.ySize; i++){
		openList.insert(std::make_pair(0.0, std::make_pair(i,0)));
	}

	// set destination flag to false

	bool foundDest = false;

	// begin loop to find path. If openList is empty, terminate the loop

	while(!openList.empty()){
		// First step is to pop the fist entry on the list
		OpenListInfo pop = *openList.begin();

		// remove from open list
		openList.erase(openList.begin());

		// Add to the closed list
		int row = pop.second.first; // first argument of the second pair
		int col = pop.second.second; // second argument of second pair
		closedList[row*info.xSize + col] = true;

		/*
			Now we need to generate all 4 successors from the popped cell.
			The successors are north, south, east, and west.
			
			North index = i - 1, j
			South index = i + 1, j
			East index =  i    , j + 1
			West index =  i    , j - 1
		*/
		float gNew, hNew, fNew;

		// Evaluate North
		
		int tempCol = col;
		int tempRow = row;

		// adjust North for periodic boundary condition

		if(row == 0){
			tempRow = info.ySize - 1;
		} else{
			tempRow = row - 1;
		}

		// check if we reached destination, which is the entire right boundary
		if(tempCol == info.ySize - 1 && GRID[tempRow*info.xSize + tempCol] != 1){
			nodeInfo[tempRow][tempCol].parentRow = row;
			nodeInfo[tempRow][tempCol].parentCol = col;
			if(info.verbose == true){
				printf("The destination cell was found.\n");
			}
			// Call trace path function
			foundDest = true;
			return foundDest;
		} else if(closedList[tempRow*info.xSize + tempCol] == false && GRID[tempRow*info.xSize + tempCol] == 0) // check if successor is not on closed list and not a solid wall
		{
			gNew = nodeInfo[row][col].g + 1.0;	// cost from moving from last cell to this cell
			hNew = (info.xSize - 1) - tempCol; // Since entire right boundary is the distance, h is just a count of the number of columns from the right.	
			fNew = gNew + hNew;					// total cost is just h+g
			// Check if on open list. If yes, update f,g, and h accordingly.
			// If not, add it to open list.
			if(nodeInfo[tempRow][tempCol].f == FLT_MAX || nodeInfo[tempRow][tempCol].f > fNew){
				openList.insert(std::make_pair(fNew, std::make_pair(tempRow, tempCol)));
				nodeInfo[tempRow][tempCol].f = fNew;
				nodeInfo[tempRow][tempCol].g = gNew;
				nodeInfo[tempRow][tempCol].h = hNew;
				nodeInfo[tempRow][tempCol].parentRow = row;
				nodeInfo[tempRow][tempCol].parentCol = col;
			}
		}
			

		// Evaluate South

		tempCol = col;
		tempRow = row;

		// Adjust for periodic BC

		if(row == info.ySize - 1){
			tempRow = 0;
		} else{
			tempRow = row + 1;
		}

		// check if we reached destination, which is the entire right boundary
		if(tempCol == info.ySize - 1 && GRID[tempRow*info.xSize + tempCol] != 1){
			nodeInfo[tempRow][tempCol].parentRow = row;
			nodeInfo[tempRow][tempCol].parentCol = col;
			if(info.verbose == true){
				printf("The destination cell was found.\n");
			}
			// Call trace path function
			foundDest = true;
			return foundDest;
		} else if(closedList[tempRow*info.xSize + tempCol] == false && GRID[tempRow*info.xSize + tempCol] == 0) // check if successor is not on closed list and not a solid wall
		{
			gNew = nodeInfo[row][col].g + 1.0;	// cost from moving from last cell to this cell
			hNew = (info.xSize - 1) - tempCol; // Since entire right boundary is the distance, h is just a count of the number of columns from the right.	
			fNew = gNew + hNew;					// total cost is just h+g
			// Check if on open list. If yes, update f,g, and h accordingly.
			// If not, add it to open list.
			if(nodeInfo[tempRow][tempCol].f == FLT_MAX || nodeInfo[tempRow][tempCol].f > fNew){
				openList.insert(std::make_pair(fNew, std::make_pair(tempRow, tempCol)));
				nodeInfo[tempRow][tempCol].f = fNew;
				nodeInfo[tempRow][tempCol].g = gNew;
				nodeInfo[tempRow][tempCol].h = hNew;
				nodeInfo[tempRow][tempCol].parentRow = row;
				nodeInfo[tempRow][tempCol].parentCol = col;
			}
		}

		// Evaluate East (if it exists)

		if(col != info.xSize - 1){
			tempRow = row;
			tempCol = col + 1;

			// check if we reached destination, which is the entire right boundary
			if(tempCol == info.ySize - 1 && GRID[tempRow*info.xSize + tempCol] != 1){
				nodeInfo[tempRow][tempCol].parentRow = row;
				nodeInfo[tempRow][tempCol].parentCol = col;
				if(info.verbose == true){
					printf("The destination cell was found.\n");
				}
				// Call trace path function
				foundDest = true;
				return foundDest;
			} else if(closedList[tempRow*info.xSize + tempCol] == false && GRID[tempRow*info.xSize + tempCol] == 0) // check if successor is not on closed list and not a solid wall
			{
				gNew = nodeInfo[row][col].g + 1.0;	// cost from moving from last cell to this cell
				hNew = (info.xSize - 1) - tempCol; // Since entire right boundary is the distance, h is just a count of the number of columns from the right.	
				fNew = gNew + hNew;					// total cost is just h+g
				// Check if on open list. If yes, update f,g, and h accordingly.
				// If not, add it to open list.
				if(nodeInfo[tempRow][tempCol].f == FLT_MAX || nodeInfo[tempRow][tempCol].f > fNew){
					openList.insert(std::make_pair(fNew, std::make_pair(tempRow, tempCol)));
					nodeInfo[tempRow][tempCol].f = fNew;
					nodeInfo[tempRow][tempCol].g = gNew;
					nodeInfo[tempRow][tempCol].h = hNew;
					nodeInfo[tempRow][tempCol].parentRow = row;
					nodeInfo[tempRow][tempCol].parentCol = col;
				}
			}
		}

		// Evaluate West

		if(col != 0){
			tempRow = row;
			tempCol = col;

			// check if we reached destination, which is the entire right boundary
			if(tempCol == info.ySize - 1 && GRID[tempRow*info.xSize + tempCol] != 1){
				nodeInfo[tempRow][tempCol].parentRow = row;
				nodeInfo[tempRow][tempCol].parentCol = col;
				if(info.verbose == true){
					printf("The destination cell was found.\n");
				}
				// Call trace path function
				foundDest = true;
				return foundDest;
			} else if(closedList[tempRow*info.xSize + tempCol] == false && GRID[tempRow*info.xSize + tempCol] == 0) // check if successor is not on closed list and not a solid wall
			{
				gNew = nodeInfo[row][col].g + 1.0;	// cost from moving from last cell to this cell
				hNew = (info.xSize - 1) - tempCol; // Since entire right boundary is the distance, h is just a count of the number of columns from the right.	
				fNew = gNew + hNew;					// total cost is just h+g
				// Check if on open list. If yes, update f,g, and h accordingly.
				// If not, add it to open list.
				if(nodeInfo[tempRow][tempCol].f == FLT_MAX || nodeInfo[tempRow][tempCol].f > fNew){
					openList.insert(std::make_pair(fNew, std::make_pair(tempRow, tempCol)));
					nodeInfo[tempRow][tempCol].f = fNew;
					nodeInfo[tempRow][tempCol].g = gNew;
					nodeInfo[tempRow][tempCol].h = hNew;
					nodeInfo[tempRow][tempCol].parentRow = row;
					nodeInfo[tempRow][tempCol].parentCol = col;
				}
			}
		}
	}
	if(info.verbose == true){
		printf("Failed to find a path.\n");
	}
	return foundDest;
}


int explicitMomentum(float *uExp, float *vExp, float *u, float *v, float *uCoeff, float *vCoeff,
	options* o, simulationInfo* info)
{

}