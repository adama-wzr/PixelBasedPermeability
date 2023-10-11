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
	float alphaRelax;
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
	printf("Relaxation Factor = %f\n", opts->alphaRelax);
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

	 	}else if(strcmp(tempC, "RelaxFactor:") == 0){
	 		opts->alphaRelax = tempD;

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


float findMax(float one, float two, float three){
	/*
		findMax Function

		Inputs:
			- three floats
		Outputs:
			- returns the maximum of the three floats.
	*/

	float max;
	float target[3] = {one, two, three};
	for(int i = 0; i<3; i++){
		if(i == 0){
			max = target[i];
		} else if(target[i] > max){
			max = target[i];
		}
	}

	return max;
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


int explicitMomentum(unsigned int *Grid, float *uExp, float *vExp, float *u, float *v, float *uCoeff, float *vCoeff,
	options* o, simulationInfo* info)
{

	/*
		Explicit Momentum:

		Inputs:
			- unsigned int *Grid: pointer to array holding boundary information
			- float *uExp: pointer to array to store the explicit component of u-velocity
			- float *vExp: pointer to array to store the explicit component of v-velocity
			- float *u: pointer to array storing u-velocity (from last iterative step)
			- float *v: pointer to array storing v-velocity (from last iterative step)
			- float *uCoeff: pointer to array to store the central coefficient of every u-velocity
			- float *vCoeff: pointer to array to store the central coefficient of every v-velocity
			- options* o: pointer to options data structure
			- simulationInfo* info: pointer to simulationInfo data-structure
		Outpus:
			None.

		Function will explicitly solve and modify the arrays uExp and vExp accordingly. Since this is
		explicit, the discretization and solution are done in the same loop. 

	*/

	// Local variables for dx and dy and Area of each cell
	float dx, dy;
	dx = info->dx;
	dy = info->dy;
	float alpha = o->alphaRelax;
	float A = dx*dy;
	// Useful variables in the solution
	float fwU, fwV, feU, feV, fnU, fnV, fsU, fsV;
	float sourceLHS, sourceRHS;
	float DeltaF;
	float awU, awV, aeU, aeV, anU, anV, asU, asV;
	float density = o->density;
	float viscosity = o->viscosity;

	float temp[3];

	// Useful indexing variables

	int nColsU = info->numCellsX + 1;
	int nColsV = info->numCellsY;

	int nRowsV = info->numCellsY+1;
	int nRowsU = info->numCellsY;

	int uRow, uCol, vRow, vCol;

	// Define d: only possible because visc = constant and grid is uniform

	float d = viscosity*A/dx;

	// Write discretization
	for(int i = 0; i<info->numCellsX+1; i++){
		for(int j = 0; j<info->numCellsY; j++){
			uRow = j;
			uCol = i;

			vRow = i;
			vCol = j;

			// Explicit U coefficients

			// fw and fe don't depend on corners
			if (uCol == 0){
				fwU = 1/2*density*A*u[uRow*nColsU + uCol];
				feU = 1/2*density*A*(u[uRow*nColsU + uCol] + u[uRow*nColsU + uCol + 1]);
			} else if(uCol == info->numCellsX){
				fwU = 1/2*density*A*(u[uRow*nColsU + uCol] + u[uRow*nColsU + uCol - 1]);
				feU = 1/2*density*A*u[uRow*nColsU + uCol];
			} else{
				fwU = 1/2*density*A*(u[uRow*nColsU + uCol] + u[uRow*nColsU + uCol - 1]);
				feU = 1/2*density*A*(u[uRow*nColsU + uCol] + u[uRow*nColsU + uCol + 1]);
			}

			// fs and fn depend on boundaries a lot more
			if(uCol == 0 && uRow == 0){
				// top left corner
				fnU = 1/2*density*A*v[(nRowsV - 1)*nColsV + uCol];
				fsU = 1/2*density*A*v[(uRow + 1)*nColsV + uCol];
			} else if(uCol == 0 && uRow == nRowsU - 1){
				// bottom left
				fsU = 1/2*density*A*v[0];
				fnU = 1/2*density*A*v[(uRow+1)*nColsV + uCol];
			} else if(uCol == 0){
				// Anywhere in left boundary
				fnU = 1/2*density*A*v[uRow*nColsV + uCol];
				fsU = 1/2*density*A*v[(uRow + 1)*nColsV + uCol];
			} else if(uCol == nColsU - 1 && uRow == 0){
				// top right corner
				fnU = 1/2*density*A*v[(nRowsV - 1)*nColsV + uCol - 1];
				fsU = 1/2*density*A*v[(uRow + 1)*nColsV + uCol - 1];
			}else if(uCol == nColsU - 1 && uRow == nRowsU - 1){
				// bottom right corner
				fnU = 1/2*density*A*v[(uRow)*nColsV + uCol - 1];
				fsU = 1/2*density*A*v[(0)*nColsV + uCol - 1];
			} else if(uCol == nColsU - 1){
				// right boundary
				fnU = 1/2*density*A*v[(uRow)*nColsV + uCol - 1];
				fsU = 1/2*density*A*v[(uRow + 1)*nColsV + uCol - 1];
			} else if(uRow == 0){
				// top boundary
				fnU = 1/2*density*A*(v[(nRowsV - 1)*nColsV + uCol] + v[(nRowsV - 1)*nColsV + uCol - 1]);
				fsU = 1/2*density*A*(v[(uRow + 1)*nColsV + uCol] + v[(uRow + 1)*nColsV + uCol - 1]);
			} else if(uRow == nRowsU - 1){
				// bottom boundary
				fnU = 1/2*density*A*(v[uRow*nColsV + uCol] + v[uRow*nColsV + uCol - 1]);
				fsU = 1/2*density*A*(v[uCol] + v[uCol - 1]);
			} else{
				// not a boundary
				fsU = 1/2*density*A*(v[(uRow + 1)*nColsV + uCol] + v[(uRow + 1)*nColsV + uCol - 1]);
				fnU = 1/2*density*A*(v[uRow*nColsV + uCol] + v[uRow*nColsV + uCol - 1]);
			}
			

			// Explicit V-coefficients

			if (vRow == 0){
				// top
				fwV = 1/2*density*A*(u[vRow*nColsU + vCol] + u[(nRowsU - 1)*nColsU + vCol]);
				feV = 1/2*density*A*(u[vRow*nColsU + vCol + 1] + u[(nRowsU - 1)*nColsU + vCol + 1]);

				fsV = 1/2*density*A*(v[vRow*nColsV + vCol] + v[(vRow + 1)*nColsV + vCol]);
				fnV = 1/2*density*A*(v[vRow*nColsV + vCol] + v[(nRowsV - 1)*nColsV + vCol]);
			} else if(vRow == nRowsV - 1){
				// bottom
				fwV = 1/2*density*A*(u[(vRow - 1)*nColsU + vCol] + u[(0)*nColsU + vCol]);
				feV = 1/2*density*A*(u[(vRow - 1)*nColsU + vCol + 1] + u[(0)*nColsU + vCol + 1]);

				fsV = 1/2*density*A*(v[vRow*nColsV + vCol] + v[vCol]);
				fnV = 1/2*density*A*(v[vRow*nColsV + vCol] + v[(vRow - 1)*nColsV + vCol]);
			} else{
				// not a boundary
				fwV = 1/2*density*A*(u[vRow*nColsU + vCol] + u[(vRow - 1)*nColsU + vCol]);
				feV = 1/2*density*A*(u[vRow*nColsU + vCol + 1] + u[(vRow - 1)*nColsU + vCol + 1]);

				fsV = 1/2*density*A*(v[vRow*nColsV + vCol] + v[(vRow + 1)*nColsV + vCol]);
				fnV = 1/2*density*A*(v[vRow*nColsV + vCol] + v[(vRow - 1)*nColsV + vCol]);
			}

			// Check if U is in a solid interface. If not gather coefficients and calculate explicit component

			if(uCol == 0 && Grid[uRow*nColsV + uCol] == 1){
				uExp[uRow*nColsU + uCol] = 0;
				uCoeff[uRow*nColsU + uCol] = 1;
			} else if(uCol == nColsU - 1 && Grid[uRow*nColsV + uCol - 1] == 1){
				uExp[uRow*nColsU + uCol] = 0;
				uCoeff[uRow*nColsU + uCol] = 1;
			} else if(Grid[uRow*nColsV + uCol] == 1 || Grid[uRow*nColsV + uCol - 1] == 1){
				uExp[uRow*nColsU + uCol] = 0;
				uCoeff[uRow*nColsU + uCol] = 1;
			} else{
				// Hybrid Discretization
				awU = findMax(fwU, d + fwU/2, 0.0f);
				aeU = findMax(-feU, d - feU/2, 0.0f);
				asU = findMax(fsU, d + fsU/2, 0.0f);
				anU = findMax(-fnU, d - fnU/2, 0.0f);

				DeltaF = feU - fwU + fnU - fsU;

				// Store central coefficient

				// Should we select these coefficient better according to boundaries?

				uCoeff[uRow*nColsU + uCol] = awU + aeU + asU + anU + DeltaF;

				// Initalize uExplicit as a component of the previous iterative step

				uExp[uRow*nColsU + uCol] = (1 - o->alphaRelax)*u[uRow*nColsU + uCol];

				// Increment uExp based on neighborhood

				if(uCol != 0){
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(awU*u[uRow*nColsU + uCol - 1]);
				}

				if(uCol != nColsU - 1){
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(aeU*u[uRow*nColsU + uCol + 1]);
				}

				if(uRow != 0){
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(anU*u[(uRow - 1)*nColsU + uCol]);
				}

				if(uRow != nRowsU - 1){
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(awU*u[(uRow + 1)*nColsU + uCol]);
				}
			}

			// Check if V is in a solid interface. If not gather coefficients and calculate explicit component

			if(vRow == 0 && Grid[vRow*nColsV + vCol] == 1){
				vExp[vRow*nColsV + vCol] = 0;
				vCoeff[vRow*nColsV + vCol] = 1;
			} else if(vRow == nRowsV - 1 && Grid[(vRow - 1)*nColsV + vCol] == 1){
				vExp[vRow*nColsV + vCol] = 0;
				vCoeff[vRow*nColsV + vCol] = 1;
			} else if(Grid[vRow*nColsV + vCol] == 1 || Grid[(vRow - 1)*nColsV + vCol] == 1){
				vExp[vRow*nColsV + vCol] = 0;
				vCoeff[vRow*nColsV + vCol] = 1;
			} else{
				// Hybrid Discretization
				awV = findMax(fwV, d + fwV/2, 0.0f);
				aeV = findMax(-feV, d - feV/2, 0.0f);
				asV = findMax(fsV, d + fsV/2, 0.0f);
				anV = findMax(-fnV, d - fnV/2, 0.0f);

				DeltaF = feV - fwV + fnV - fsV;

				// Store central coefficient

				// Should we select these coefficient better according to boundaries?

				vCoeff[vRow*nColsV + vCol] = awV + aeV + asV + anV + DeltaF;

				// Initalize vExplicit as a component of the previous iterative step

				vExp[vRow*nColsV + vCol] = (1 - o->alphaRelax)*v[vRow*nColsV + vCol];

				// Increment vExp based on neighborhood

				if(vCol != 0){
					vExp[vRow*nColsV + vCol] = o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(awV*v[vRow*nColsV + vCol - 1]);
				}

				if(vCol != nColsV-1){
					vExp[vRow*nColsV + vCol] = o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(aeV*v[vRow*nColsV + vCol + 1]);
				}

				if(vRow != 0){
					vExp[vRow*nColsV + vCol] = o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(anV*v[(vRow - 1)*nColsV + vCol]);
				}

				if(vRow != nRowsV - 1){
					vExp[vRow*nColsV + vCol] = o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(asV*v[(vRow + 1)*nColsV + vCol]);
				}
			} // end of if statement

		}
	} // end of for loops

	return 0;
}

int implicitPressure(unsigned int *Grid, float *uExp, float *vExp, float *uCoeff, float *vCoeff, float *Pressure,
	options* o, simulationInfo* info)
{
	/*
	Implicit Pressure:

	Inputs:
		- unsigned int *Grid: pointer to array holding boundary information
		- float *uExp: pointer to array to store the explicit component of u-velocity
		- float *vExp: pointer to array to store the explicit component of v-velocity
		- float *uCoeff: pointer to array to store the central coefficient of every u-velocity
		- float *vCoeff: pointer to array to store the central coefficient of every v-velocity
		- float *Pressure: array containing the pressure at each step
		- options* o: pointer to options data structure
		- simulationInfo* info: pointer to simulationInfo data-structure
	Outpus:
		None.

	Function will assemble a coefficient matrix and a RHS to implicitly solve for Pressure. The Pressure array is
	modified with the final answer at this iterative step. 

	*/

	float dx, dy;
	dx = info->dx;
	dy = info->dy;
	float alpha = o->alphaRelax;
	float Area = dx*dy;
	float density = o->density;
	float viscosity = o->viscosity;

	/* Create the "A" matrix in Ax = b system
	
	Indexing is as follows:
	
	P = 0
	W = 1
	E = 2
	S = 3
	N = 4
	*/

	float *A = (float *)malloc(sizeof(float)*info->numCellsX*info->numCellsY*5);
	float *RHS = (float *)malloc(sizeof(float)*info->numCellsX*info->numCellsY);

	// Initialize A to all zero's to avoid NaN's
	memset(A, 0, sizeof(A));

	// Indexing variables

	int nColsU = info->numCellsX + 1;
	int nColsV = info->numCellsX;
	int nColsP = info->numCellsX;
	int index = 0;

	// Neighborhood parameters

	float ae, aw, as, an;
	float aE, aW, aS, aN, aP;
	float bp;

	for(int row = 0; row<info->numCellsY; row++){
		for(int col = 0; col<info->numCellsX; col++){
			index = row*nColsP + col;

			// Initialize to 0's to avoid NaN's
			aP = 0;
			aN = 0;
			aS = 0;
			aW = 0;
			aE = 0;
			bp = 0;
			RHS[index] = 0;

			// Check if P[index] is solid. If so, equation is P = 0 (central coeff = 1, RHS = 0)
			if(Grid[index] == 1){
				A[index*5 + 0] = 1;
			}else
			{ 
				// If central point is not solid, then we proceed normally
				if (row!=0)
				{
					// This means a North exists
					if(Grid[index - nColsP] == 1)
					{
						// North is a solid
						aN = 0;
					} else
					{
						// North is not a solid
						aN = alpha*density*Area*Area/vCoeff[row*nColsV + col];
						aP += aN;
						bp += density*vExp[row*nColsV + col]*Area;
					}
				} else
				{
					// row == 0, apply periodic BC
					if(Grid[(info->numCellsY - 1)*nColsP + col] == 1)
					{
						// Periodic North is solid
						aN = 0;
					} else
					{
						// Periodic North is not a solid, proceed normally
						aN = alpha*density*Area*Area/vCoeff[(info->numCellsY)*nColsV + col];
						aP += aN;
						bp += density*vExp[(info->numCellsY)*nColsV + col]*Area;
					}
				}

				// Check South
				if(row != info->numCellsY - 1)
				{
					// this means south exists
					if (Grid[index + nColsP] == 1)
					{
						// South is a solid
						aS = 0;
					}else
					{
						// South exists and is not a solid
						aS = alpha*density*Area*Area/vCoeff[(row + 1)*nColsV + col];
						aP += aS;
						bp += -density*vExp[(row + 1)*nColsV + col]*Area;
					}
				} else
				{
					// Apply periodic BC
					if(Grid[col] == 1)
					{
						// This means it is solid
						aS = 0;
					}else
					{
						// Not solid
						aS = alpha*density*Area*Area/vCoeff[col];
						aP += aS;
						bp += -density*vExp[col]*Area;
					}
				}

				// aE and aW always exist, as long as they are not solid
				if(col == 0)
				{
					// We are at the LHS boundary and we know our grid is not solid
					aW = alpha*density*Area*Area/uCoeff[row*nColsU + col];
					bp += -density*uExp[row*nColsU + col]*Area;
					RHS[index] += -aW*o->PL;
					aP += aW;
					// Check that east grid is not solid. If it is set coeff to zero
					if(Grid[index + 1] == 1)
					{
						aE = 0;
					}else
					{
						aE = alpha*density*Area*Area/uCoeff[row*nColsU + col + 1];
						bp += density*uExp[row*nColsU + col + 1]*Area;
						aP += aE;
					}
				} else if(col == info->numCellsX - 1)
				{
					// We are at the RHS boundary and we know our grid is not solid
					aE = alpha*density*Area*Area/uCoeff[row*nColsU + col + 1];
					bp += density*uExp[row*nColsU + col + 1]*Area;
					RHS[index] += -aE*o->PR;
					aP += aE;
					// Check that west grid is not solid. If it is set coeff to zero
					if(Grid[index - 1] == 1)
					{
						aW = 0;
					} else
					{
						aW = alpha*density*Area*Area/uCoeff[row*nColsU + col];
						bp += -density*uExp[row*nColsU + col]*Area;
						aP += aW;
					}
				} else
				{
					// we are not in any boundary
					// Check if West is solid
					if(Grid[index - 1] == 1)
					{
						aW = 0;
					}else
					{
						aW = alpha*density*Area*Area/uCoeff[row*nColsU + col];
						bp += -density*uExp[row*nColsU + col]*Area;
						aP += aW;
					}
					// Do thhe same for East
					if(Grid[index + 1] == 1)
					{
						aE = 0;
					} else
					{
						aE = alpha*density*Area*Area/uCoeff[row*nColsU + col + 1];
						bp += density*uExp[row*nColsU + col + 1]*Area;
						aP += aE;
					}
				}
				// Gather Coefficients into the linear system of eqs.
				RHS[index] += bp;

				A[index*5 + 0] = -aP;
				A[index*5 + 1] = aW;
				A[index*5 + 2] = aE;
				A[index*5 + 3] = aS;
				A[index*5 + 4] = aN;
			}
		}
	}

	// Solve

	// Return

	return 0;

}