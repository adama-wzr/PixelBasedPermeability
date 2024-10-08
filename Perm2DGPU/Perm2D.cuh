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
	int nCores;
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
	float Flowrate;
	float ResidualX;
	float ResidualY;
	float ResidualP;
}simulationInfo;

// define pair for coords

typedef std::pair<int, int> coordPair;


__global__ void updateX_V1(float* A, float* x, float* b, float* xNew, int numCellsX, int numCellsY)
{
	// Figure out the index
	unsigned int myIdx = blockIdx.x * blockDim.x + threadIdx.x;
	// break index into row and col
	int myRow = myIdx/numCellsX;
	int myCol = myIdx % numCellsX;

	int nElements = numCellsX*numCellsY;

	// do the Jacobi Iteration with periodic BC

	if (myIdx < nElements){
		float sigma = 0;
		for(int j = 1; j<5; j++){
			if(A[myIdx*5 + j] !=0){
				if(j == 1){
					sigma += A[myIdx*5 + j]*x[myIdx - 1];
				} else if(j == 2){
					sigma += A[myIdx*5 + j]*x[myIdx + 1];
				} else if(j == 3){
					if(myRow == numCellsY - 1){
						sigma += A[myIdx*5 + j]*x[myCol];
					}else{
						sigma += A[myIdx*5 + j]*x[myIdx + numCellsX];
					}
				} else if(j == 4){
					if(myRow == 0){
						sigma += A[myIdx*5 + j]*x[(numCellsY - 1)*numCellsX + myCol];
					}else{
						sigma += A[myIdx*5 + j]*x[myIdx - numCellsX];
					}
				}
			}
		}
		xNew[myIdx] = 1/A[myIdx*5 + 0] * (b[myIdx] - sigma);
	}
}


__global__ void updateX_SOR(float* A, float* x, float* b, float* xNew, int numCellsX, int numCellsY)
{
	// Figure out the index
	unsigned int myIdx = blockIdx.x * blockDim.x + threadIdx.x;
	// break index into row and col
	int myRow = myIdx/numCellsX;
	int myCol = myIdx % numCellsX;
	float w = 2.0/3.0;

	int nElements = numCellsX*numCellsY;

	// do the Jacobi Iteration with periodic BC

	if (myIdx < nElements){
		float sigma = 0;
		for(int j = 1; j<5; j++){
			if(A[myIdx*5 + j] !=0){
				if(j == 1){
					sigma += A[myIdx*5 + j]*x[myIdx - 1];
				} else if(j == 2){
					sigma += A[myIdx*5 + j]*x[myIdx + 1];
				} else if(j == 3){
					if(myRow == numCellsY - 1){
						sigma += A[myIdx*5 + j]*x[myCol];
					}else{
						sigma += A[myIdx*5 + j]*x[myIdx + numCellsX];
					}
				} else if(j == 4){
					if(myRow == 0){
						sigma += A[myIdx*5 + j]*x[(numCellsY - 1)*numCellsX + myCol];
					}else{
						sigma += A[myIdx*5 + j]*x[myIdx - numCellsX];
					}
				}
			}
		}
		xNew[myIdx] = (1.0-w)*x[myIdx] + w/A[myIdx*5 + 0] * (b[myIdx] - sigma);
	}
}



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
	printf("Dinamic Viscosity = %.6f kg/(m*s) or Pa*s\n", opts->viscosity);
	printf("Mesh Refinement = %d\n", opts->MeshAmp);
	printf("Relaxation Factor = %f\n", opts->alphaRelax);
	printf("Maximum Iterations Solver = %ld\n", opts->MaxIterSolver);
	printf("Solver Convergence = %.10f\n", opts->ConvergenceSolver);
	printf("Maximum Iterations Global = %ld\n", opts->MaxIterGlobal);
	printf("RMS Convergence = %.10f\n", opts->ConvergenceRMS);
	printf("Name of input image: %s\n", opts->inputFilename);
	printf("Name of output file: %s\n", opts->outputFilename);
	printf("Number of Cores = %d\n", opts->nCores);
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

	 	}else if(strcmp(tempC, "nCores:") == 0){
	 		opts->nCores = (int)tempD;
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


int initializeGPU(float **d_x_vec, float **d_temp_x_vec, float **d_RHS, float **d_Coeff, simulationInfo* mesh){

	// Set device, when cudaStatus is called give status of assigned device.
	// This is important to know if we are running out of GPU space
	cudaError_t cudaStatus = cudaSetDevice(0);

	// Start by allocating space in GPU memory

	if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		getchar();
        return 0;
    }

    cudaStatus = cudaMalloc((void**)&(*d_x_vec), mesh->nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		getchar();
        return 0;
    }

    cudaStatus = cudaMalloc((void**)&(*d_temp_x_vec), mesh->nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		getchar();
        return 0;
    }

    cudaStatus = cudaMalloc((void**)&(*d_RHS), mesh->nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		getchar();
        return 0;
    }

    cudaStatus = cudaMalloc((void**)&(*d_Coeff), mesh->nElements*sizeof(float)*5);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		getchar();
        return 0;
    }

    // Set GPU buffers (initializing matrices to 0)

     // Memset GPU buffers
    cudaStatus = cudaMemset((*d_x_vec),0, mesh->nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemset failed!");
		getchar();
        return 0;
    }

	// Memset GPU buffers
    cudaStatus = cudaMemset((*d_temp_x_vec),0, mesh->nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemset failed!");
		getchar();
        return 0;
    }

     // Memset GPU buffers
    cudaStatus = cudaMemset((*d_RHS),0, mesh->nElements*sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemset failed!");
		getchar();
        return 0;
    }

	// Memset GPU buffers
    cudaStatus = cudaMemset((*d_Coeff),0, 5*mesh->nElements*sizeof(float));		// coefficient matrix has the 5 main diagonals for all elements
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemset failed!");
		getchar();
        return 0;
    }

    return 1;
}


void unInitializeGPU(float **d_x_vec, float **d_temp_x_vec, float **d_RHS, float **d_Coeff)
{
	cudaError_t cudaStatus;

	if((*d_x_vec)!=NULL)
    cudaStatus = cudaFree((*d_x_vec));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaFree failed!");
        return;
    }

	if((*d_temp_x_vec)!=NULL)
    cudaStatus = cudaFree((*d_temp_x_vec));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaFree failed!");
        return;
    }

	if((*d_Coeff)!=NULL)
    cudaStatus = cudaFree((*d_Coeff));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaFree failed!");
        return;
    }

	if((*d_RHS)!=NULL)
    cudaStatus = cudaFree((*d_RHS));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaFree failed!");
        return;
    }    

	cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
		getchar();
        return;
    }
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


int printPUVmaps(float* Pressure, float* u, float* v, options *o, simulationInfo *info){
	/*
		printPUVmaps
		Inputs:
			- float* Pressure: pointer to array containing Pressure data
			- float* u: pointer to array containing the x-component of velocity
			- float* v: pointer to array containing the y-component of velocity
			- options* o: pointer to datastructure containing user entered options
			- simulationInfo *info: pointer to datastructure containing simulation info

		Output:
			- None

		Function prints maps to a .csv file (user entered filename).
	*/

	FILE *MAP;
	MAP = fopen(o->outputFilename, "w");

	fprintf(MAP, "P,U,V,x,y\n");

	int nColsU = info->numCellsX + 1;
	int nColsV = info->numCellsX;
	int nColsP = info->numCellsX;

	for(int i = 0; i<info->numCellsY; i++){
		for(int j = 0; j<info->numCellsX; j++){
			float uc = (u[i*nColsU + j] + u[i*nColsU + j + 1])/2;
			float vc = (v[(i + 1)*nColsV + j] + v[i*nColsV + j])/2;
			fprintf(MAP,"%f,%f,%f,%d,%d\n", Pressure[i*nColsP + j], uc, vc, j, i);
		}
	}

	fclose(MAP);

	return 0;
}


float ResMap(float *U, float *V, options *o, simulationInfo *info){
	/*
		Funcion ResMap:

		Inputs:
			- float *U: pointer to array with x component of velocities
			- float *V: pointer to array with v components of velocities
			- options *o: pointer to array with user-entered options
			- simulationInfo *info: pointer to array with simulation domain information

		Outputs:
			- Residual

		Function will calculate residual convergence in the continuity equation, and each iterative step is normalized by 
		the absolute value of the RHS summed over all cells for that iterative step.
	*/

	// Domain variables

	float dx, dy;
	dx = info->dx;
	dy = info->dy;
	float Area = dx*dy;

	// Useful indexing variables

	int nColsU = info->numCellsX + 1;
	int nColsV = info->numCellsY;

	// Properties

	float density = o->density;

	// Iterate through the entire domain

	float R = 0; 	// variable used to store the residual

	float max = 0;

	float cellCont = 0;

	FILE *MAP;
	MAP = fopen("ResMap.csv", "w");

	fprintf(MAP, "R,x,y\n");

	for(int row = 0; row<info->numCellsY; row++){
		for(int col = 0; col<info->numCellsX; col++){
			cellCont = density*Area*(fabs(U[row*nColsU + col] - U[row*nColsU + col + 1] + V[(row + 1)*nColsV + col] - V[row*nColsV + col]));
			fprintf(MAP, "%1.12f,%d,%d\n", cellCont, col, row);
			R += cellCont;
			if(cellCont > max){
				max = cellCont;
			}
		}
	}

	printf("Max Cell Continuity = %1.9f\n", max);
	fclose(MAP);
	// R = R/(info->numCellsY*info->numCellsX);

	return max;
}


float ResidualContinuity(float *U, float *V, options *o, simulationInfo *info){
	/*
		Funcion ResidualContinuity:

		Inputs:
			- float *U: pointer to array with x component of velocities
			- float *V: pointer to array with v components of velocities
			- options *o: pointer to array with user-entered options
			- simulationInfo *info: pointer to array with simulation domain information

		Outputs:
			- Residual

		Function will calculate residual convergence in the continuity equation, and each iterative step is normalized by 
		the absolute value of the RHS summed over all cells for that iterative step.
	*/

	// Domain variables

	float dx, dy;
	dx = info->dx;
	dy = info->dy;
	float Area = dx*dy;

	// Useful indexing variables

	int nColsU = info->numCellsX + 1;
	int nColsV = info->numCellsY;

	// Properties

	float density = o->density;

	// Iterate through the entire domain

	float R = 0; 	// variable used to store the residual

	float max = 0;

	float cellCont = 0;

	for(int row = 0; row<info->numCellsY; row++){
		for(int col = 1; col<info->numCellsX - 1; col++){
			cellCont = density*Area*(fabs(U[row*nColsU + col] - U[row*nColsU + col + 1] + V[(row + 1)*nColsV + col] - V[row*nColsV + col]));
			R += cellCont;
			if(cellCont > max){
				max = cellCont;
			}
		}
	}

	printf("Max Cell Continuity = %e\n", max);
	// R = R/(info->numCellsY*info->numCellsX);

	return max;
}


int PermCalc(float *U, options *o, simulationInfo *info){
	/*
		Funcion PermCalc:

		Inputs:
			- float *U: pointer to array with x component of velocities
			- options *o: pointer to array with user-entered options
			- simulationInfo *info: pointer to array with simulation domain information

		Outputs:
			- None

		Function will calculate the normalized permeability at each step.
	*/

	float dx, dy;
	dx = info->dx;
	dy = info->dy;
	float Area = dx*dy;	// this is the cross-sectional area of each cell
	float viscosity = o->viscosity;

	int nRowsU = info->numCellsY;
	int nColsU = info->numCellsX + 1;

	float *FlowRate = (float *)malloc(sizeof(float)*nColsU);

	for(int k = 0; k<nColsU; k++){
		FlowRate[k] = 0;
		for(int row = 0; row<nRowsU; row++){
			FlowRate[k] += U[row*nColsU + k]*Area;	
		}
	}

	float Qavg = 0;

	for(int k = 0; k<nColsU; k++){
		Qavg += FlowRate[k];
	}

	Qavg = Qavg/nColsU;

	free(FlowRate);

	info->Flowrate = Qavg;

	// info->Perm = Qavg/(o->DomainHeight*dx)*viscosity*o->DomainWidth/((o->PL - o->PR));
	info->Perm = Qavg/(o->DomainHeight*o->DomainWidth)*viscosity*o->DomainWidth/((o->PL - o->PR)*o->DomainHeight*dx);
	return 0;
}


float JacobiConv(float *arr, float *RHS, float *Pressure, simulationInfo *info)
{
	/*
		Function JacobiConv:

		Inputs:
			- float *arr: pointer to coefficient matrix array
			- float *sol: pointer to RHS array (Ax = b, this is b)
			- float *Pressure: pointer to x array (Ax = b, this is x)
			- simulationInfo *info: pointer to data struct with information about the mesh
		
		Outputs:
			- Returns the residual value calculated as a float

		Function will calculate a residual for assessing Jacobi convergence. The Ax = b system, in theory,
		should yield (b - Ax) = 0. In practice, the solution yields x_bar, such that (b - Ax_bar) > 0.

		This function returns r = 1/N*(sum(fabs(b - Ax_bar))) for a A being N x N and 
		x_bar and b being N in length. r is a positive scalar.
	
	*/

	float ResidualCalc = 0;

	int nRows = info->numCellsY;
	int nCols = info->numCellsX;

	for(int i = 0; i<nRows; i++){
		for(int j = 0; j<nCols; j++){
			int index = i*nCols + j;
			float sumAx = 0;

			for(int A_col = 0; A_col < 5; A_col++){
				if(arr[index*5 + A_col] != 0){
					if(A_col == 0)
					{
						sumAx += arr[index*5 + A_col]*Pressure[index];
					} else if(A_col == 1)
					{
						sumAx += arr[index*5 + A_col]*Pressure[index - 1];		
					} else if(A_col == 2)
					{
						sumAx += arr[index*5 + A_col]*Pressure[index + 1];		
					}else if(j == 3)
					{
						if(i == nRows - 1)
						{
							sumAx += arr[index*5 + A_col]*Pressure[j];
						} else
						{
							sumAx += arr[index*5 + A_col]*Pressure[index + nCols];
						}
					} else if(j == 4)
					{
						if(i == 0)
						{
							sumAx += arr[index*5 + A_col]*Pressure[(nRows - 1)*nCols + j];
						} else
						{
							sumAx += arr[index*5 + A_col]*Pressure[index - nCols];
						}
					}
				}
			}

			ResidualCalc += (1.0/(float)info->nElements)*fabs(RHS[index] - sumAx);
		}
	}

	return ResidualCalc;
}


int pJacobiCPU2D(float *arr, float *sol, float *Pressure, options *o, simulationInfo *info){

	/*
		Function pJacobiCPU2D
		
		Inputs:
			- float *arr: pointer to array containing the coefficient matrix
			- float *sol: pointer to the RHS of the discretization
			- float *Pressure: pointer to the Pressure array
			- options *o: pointer to datastructure of user entered options
			- simulationInfo *info: datastructure of simulation and domain parameters

		Outputs:
			- None

		Function will solved the linear system of equations Arr*Pressure = RHS for Pressure using
		the Jacobi Iteration method. The coefficient matrix (arr) is only 5 diagonals. Note that at the top
		and bottom of the domain, the boundary is periodic. This means an adjustment must be made
		to the coefficient matrix to make sure we are multiplying the correct Pressures.

		This algorithm also uses a simple openMP parallelization.
	*/

	/*
	
	Indexing is as follows:
	
	P = 0
	W = 1
	E = 2
	S = 3
	N = 4


	*/
	// Declare useful variables
	int iterationCount = 0;
	int iterLimit = o->MaxIterSolver;
	float tolerance = o->ConvergenceSolver;
	int nCols = info->numCellsX;
	int nRows = info->numCellsY;

	// Temporary Pressure Array

	float *tempP = (float *)malloc(sizeof(float)*nCols*nRows);

	for(int i = 0; i<nCols*nRows; i++){
		tempP[i] = Pressure[i];
	}

	// More runtime related variables
	
	float sigma = 0.0;
	float norm_diff;
	float convergence_criteria = 1;
	int index;

	while(convergence_criteria > tolerance && iterationCount < iterLimit)
	{
		#pragma omp parallel private(index, sigma)
		#pragma omp for
		for(index = 0; index<nCols*nRows; index++)
		{
			int myRow = index/nCols;
			int myCol = index % nCols;
			sigma = 0;
			for(int j = 1; j<5; j++)
			{
				if(arr[index*5 + j] != 0)
				{
					if(j == 1)
					{
						sigma += arr[index*5 + j]*tempP[index - 1];
					}else if(j == 2)
					{
						sigma += arr[index*5 + j]*tempP[index + 1];
					} else if(j == 3)
					{
						if(myRow == nRows - 1)
						{
							sigma += arr[index*5 + j]*tempP[myCol];
						} else
						{
							sigma += arr[index*5 + j]*tempP[index + nCols];
						}
					} else if(j == 4)
					{
						if(myRow == 0)
						{
							sigma += arr[index*5 + j]*tempP[(nRows - 1)*nCols + myCol];
						} else
						{
							sigma += arr[index*5 + j]*tempP[index - nCols];
						}
					}
				}
			}
			Pressure[index] = 1.0/arr[index*5 + 0]*(sol[index] - sigma);
		}

		if(iterationCount % 100000 == 0)
		{
			norm_diff = 0;
			for(index = 0; index < nCols*nRows; index++)
			{
				norm_diff += fabs((Pressure[index] - tempP[index])/(o->PL*(nCols*nRows)));
			}
			printf("Normalized Absolute Change = %f, Jacobi TOL = %f\n", norm_diff, tolerance);
		}

		convergence_criteria = norm_diff;

		for(int i = 0; i<nCols*nRows; i++){
			tempP[i] = Pressure[i];
		}

		iterationCount++;

	}

	free(tempP);
	return 0;
}


int JacobiGPU2D(float *arr, float *sol, float *Pressure, options *o, simulationInfo *info,
	float *d_x_vec, float *d_temp_x_vec, float *d_Coeff, float *d_RHS){

	/*
		Function JacobiGPU2D
		
		Inputs:
			- float *arr: pointer to array containing the coefficient matrix
			- float *sol: pointer to the RHS of the discretization
			- float *Pressure: pointer to the Pressure array
			- options *o: pointer to datastructure of user entered options
			- simulationInfo *info: datastructure of simulation and domain parameters

		Outputs:
			- None

		Function will solved the linear system of equations Arr*Pressure = RHS for Pressure using
		the Jacobi Iteration method. The coefficient matrix (arr) is only 5 diagonals. Note that at the top
		and bottom of the domain, the boundary is periodic. This means an adjustment must be made
		to the coefficient matrix to make sure we are multiplying the correct Pressures.

		This algorithm also uses a simple openMP parallelization.
	*/

	/*
	
	Indexing is as follows:
	
	P = 0
	W = 1
	E = 2
	S = 3
	N = 4


	*/
	// Declare useful variables
	int iterationCount = 0;
	int iterLimit = o->MaxIterSolver;
	float tolerance = o->ConvergenceSolver;
	int nCols = info->numCellsX;
	int nRows = info->numCellsY;
	int threads_per_block = 160;
	int numBlocks = info->nElements/threads_per_block + 1;

	const char *str = (char*) malloc(1024); // To store error string

	// Temporary Pressure Array

	float *tempP = (float *)malloc(sizeof(float)*nCols*nRows);

	for(int i = 0; i<nCols*nRows; i++){
		tempP[i] = Pressure[i];
	}

	cudaSetDevice(0);

	//Copy arrays into GPU memory

	cudaError_t cudaStatus = cudaMemcpy(d_temp_x_vec, tempP, sizeof(float) * info->nElements, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Pressure cudaMemcpy failed!");
		str = cudaGetErrorString(cudaStatus);
		fprintf(stderr, "CUDA Error!:: %s\n", str);
	}
	cudaStatus = cudaMemcpy(d_RHS, sol, sizeof(float)*info->nElements, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "d_RHS cudaMemcpy failed!");
		str = cudaGetErrorString(cudaStatus);
		fprintf(stderr, "CUDA Error!:: %s\n", str);
	}
	cudaStatus = cudaMemcpy(d_Coeff, arr, sizeof(float)*info->nElements*5, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "d_Coeff cudaMemcpy failed!");
		str = cudaGetErrorString(cudaStatus);
		fprintf(stderr, "CUDA Error!:: %s\n", str);
	}

	// More runtime related variables
	
	float norm_diff;
	float convergence_criteria = 1;
	int index;

	// Cuda event recording variables

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start, 0);

	while(convergence_criteria > tolerance && iterationCount < iterLimit)
	{
		
		// Call Kernel to Calculate new x-vector
		
		// updateX_V1<<<numBlocks, threads_per_block>>>(d_Coeff, d_temp_x_vec, d_RHS, d_x_vec, info->numCellsX, info->numCellsY);
		updateX_SOR<<<numBlocks, threads_per_block>>>(d_Coeff, d_temp_x_vec, d_RHS, d_x_vec, info->numCellsX, info->numCellsY);

		cudaDeviceSynchronize();
		
		// Convergence related material

		if(iterationCount % 1000 == 0)
		{
			cudaStatus = cudaMemcpy(Pressure, d_x_vec, sizeof(float) * info->nElements, cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Pressure cudaMemcpy failed!");
				str = cudaGetErrorString(cudaStatus);
				fprintf(stderr, "CUDA Error!:: %s\n", str);
			}
			// Copy the Pressure from last iteration from GPU to host
			cudaStatus = cudaMemcpy(tempP, d_temp_x_vec, sizeof(float) * info->nElements, cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "tempP cudaMemcpy failed!");
				str = cudaGetErrorString(cudaStatus);
				fprintf(stderr, "CUDA Error!:: %s\n", str);
			}
			norm_diff = 0;
			for(index = 0; index < nCols*nRows; index++)
			{
				norm_diff += fabs((Pressure[index] - tempP[index])/(o->PL*(nCols*nRows)));
			}

			// norm_diff = JacobiConv(arr, sol, Pressure, info);

			// Un-comment the following lines if you need more details.

			// if(iterationCount % 10000 == 0){
			// 	printf("Normalized Absolute Change = %f, Jacobi TOL = %f\n", norm_diff, tolerance);
			// }
		}

		// update temporary x-vector with new values on Device (GPU)

		cudaStatus = cudaMemcpy(d_temp_x_vec, d_x_vec, sizeof(float) * info->nElements, cudaMemcpyDeviceToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "d_temp_x_vec cudaMemcpy failed!");
			str = cudaGetErrorString(cudaStatus);
			fprintf(stderr, "CUDA Error!:: %s\n", str);
		}

		convergence_criteria = norm_diff;
		
		iterationCount++;
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);

	printf("Jacobi Convergence Time = %f, iter = %d\n", elapsedTime/1000, iterationCount);

	cudaStatus = cudaMemcpy(Pressure, d_x_vec, sizeof(float)*info->nElements, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Pressure cudaMemcpy failed!");
		str = cudaGetErrorString(cudaStatus);
		fprintf(stderr, "CUDA Error!:: %s\n", str);
	}

	free(tempP);
	return 0;
}


int FloodFill(unsigned int* Grid, simulationInfo* simInfo)
{
	/*
		Flood Fill algorithm:

		Inputs:
			- pointer to Grid: array containing the location of solids and fluids
			- simulationInfo* simInfo: pointer to struct containing the simInfo

		Outputs:
			- None

		Function will modify the Grid array. If fluid is not participating, it will
		be assigned a flag of 2.
	*/

	int* Domain = (int *)malloc(sizeof(int)*simInfo->nElements);
	int index;

	// Step 1: Initialize all solids in the domain

	for(int row = 0; row<simInfo->numCellsY; row++){
		for(int col = 0; col<simInfo->numCellsX; col++){
			index = row*simInfo->numCellsX + col;
			if(Grid[index] == 1){
				Domain[index] = 1;
			} else{
				Domain[index] = -1;
			}
		}
	}

	// Step 2: Find fluid in the two boundaries (left and right), set flag to 0, add to open list

	std::set<coordPair> cList;

	for(int row = 0; row<simInfo->numCellsY; row++){
		int indexL = row*simInfo->numCellsX;
		int indexR = row*simInfo->numCellsX + simInfo->numCellsX - 1;

		if(Domain[indexL] == -1){
			Domain[indexL] = 0;
			cList.insert(std::make_pair(row, 0));
		}

		if (Domain[indexR] == -1){
			Domain[indexR] = 0;
			cList.insert(std::make_pair(row, simInfo->numCellsX-1));
		}
	}

	while(!cList.empty())
	{
		// Pop first item in the list
		coordPair pop = *cList.begin();

		// remove from open list
		cList.erase(cList.begin());

		// Get coordinates from popped item
		int row = pop.first; // first argument of the second pair
		int col = pop.second; // second argument of second pair

		/*
			Now we need to check North, South, East, and West for more fluid:
				North: row - 1, col
				South: row + 1, col
				West: row, col - 1
				East: row, col + 1
			Details:
				- No diagonals are checked.
				- Periodic BC North and South

		*/

		int tempRow, tempCol;

		// North
		tempCol = col;

		// check periodic boundary
		if(row == 0){
			tempRow = simInfo->numCellsY - 1;
		} else{
			tempRow = row - 1;
		}

		// Update list if necessary

		if(Domain[tempRow*simInfo->numCellsX + tempCol] == -1){
			Domain[tempRow*simInfo->numCellsX + tempCol] = 0;
			cList.insert(std::make_pair(tempRow, tempCol));
		}

		// South

		tempCol = col;

		// check periodic boundary

		if(row == simInfo->numCellsY - 1){
			tempRow = 0;
		} else{
			tempRow = row + 1;
		}

		// Update list if necessary

		if(Domain[tempRow*simInfo->numCellsX + tempCol] == -1){
			Domain[tempRow*simInfo->numCellsX + tempCol] = 0;
			cList.insert(std::make_pair(tempRow, tempCol));
		}

		// West

		if(col != 0){
			tempCol = col - 1;
			tempRow = row;

			if(Domain[tempRow*simInfo->numCellsX + tempCol] == -1){
				Domain[tempRow*simInfo->numCellsX + tempCol] = 0;
				cList.insert(std::make_pair(tempRow, tempCol));
			}
		}

		// East

		if(col != simInfo->numCellsX - 1){
			tempCol = col + 1;
			tempRow = row;

			if(Domain[tempRow*simInfo->numCellsX + tempCol] == -1){
				Domain[tempRow*simInfo->numCellsX + tempCol] = 0;
				cList.insert(std::make_pair(tempRow, tempCol));
			}
		}
		// repeat everything if list is still not empty
	}

	// Every flag that is still -1 means its non-participating fluid

	for(int row = 0; row<simInfo->numCellsY; row++){
		for(int col = 0; col<simInfo->numCellsX; col++){
			index = row*simInfo->numCellsX + col;
			if(Domain[index] == -1){
				Grid[index] = 2;
			}
		}
	}

	free(Domain);

	return 0;
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
	float A = dx*dy;
	// Useful variables in the solution
	float fwU, fwV, feU, feV, fnU, fnV, fsU, fsV;
	float DeltaF;
	float awU, awV, aeU, aeV, anU, anV, asU, asV;
	float density = o->density;
	float viscosity = o->viscosity;

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
				fwU = density*A*u[uRow*nColsU + uCol];
				feU = 1.0/2*density*A*(u[uRow*nColsU + uCol] + u[uRow*nColsU + uCol + 1]);
			} else if(uCol == nColsU -1){
				fwU = 1.0/2*density*A*(u[uRow*nColsU + uCol] + u[uRow*nColsU + uCol - 1]);
				feU = density*A*u[uRow*nColsU + uCol];
			} else{
				fwU = 1.0/2*density*A*(u[uRow*nColsU + uCol] + u[uRow*nColsU + uCol - 1]);
				feU = 1.0/2*density*A*(u[uRow*nColsU + uCol] + u[uRow*nColsU + uCol + 1]);
			}

			// fs and fn depend on boundaries a lot more
			if(uCol == 0 && uRow == 0){
				// top left corner
				fnU = 1.0/2*density*A*v[(nRowsV - 1)*nColsV + uCol];
				fsU = 1.0/2*density*A*v[(uRow + 1)*nColsV + uCol];
			} else if(uCol == 0 && uRow == nRowsU - 1){
				// bottom left
				fsU = 1.0/2*density*A*v[0];
				fnU = 1.0/2*density*A*v[(uRow)*nColsV + uCol];
			} else if(uCol == 0){
				// Anywhere in left boundary
				fnU = 1.0/2*density*A*v[uRow*nColsV + uCol];
				fsU = 1.0/2*density*A*v[(uRow + 1)*nColsV + uCol];
			} else if(uCol == nColsU - 1 && uRow == 0){
				// top right corner
				fnU = 1.0/2*density*A*v[(nRowsV - 1)*nColsV + uCol - 1];
				fsU = 1.0/2*density*A*v[(uRow + 1)*nColsV + uCol - 1];
			}else if(uCol == nColsU - 1 && uRow == nRowsU - 1){
				// bottom right corner
				fnU = 1.0/2*density*A*v[(uRow)*nColsV + uCol - 1];
				fsU = 1.0/2*density*A*v[(0)*nColsV + uCol - 1];
			} else if(uCol == nColsU - 1){
				// right boundary
				fnU = 1.0/2*density*A*v[(uRow)*nColsV + uCol - 1];
				fsU = 1.0/2*density*A*v[(uRow + 1)*nColsV + uCol - 1];
			} else if(uRow == 0){
				// top boundary
				fnU = 1.0/2*density*A*(v[(nRowsV - 1)*nColsV + uCol] + v[(nRowsV - 1)*nColsV + uCol - 1]);
				fsU = 1.0/2*density*A*(v[(uRow + 1)*nColsV + uCol] + v[(uRow + 1)*nColsV + uCol - 1]);
			} else if(uRow == nRowsU - 1){
				// bottom boundary
				fnU = 1.0/2*density*A*(v[uRow*nColsV + uCol] + v[uRow*nColsV + uCol - 1]);
				fsU = 1.0/2*density*A*(v[uCol] + v[uCol - 1]);
			} else{
				// not a boundary
				fsU = 1.0/2*density*A*(v[(uRow + 1)*nColsV + uCol] + v[(uRow + 1)*nColsV + uCol - 1]);
				fnU = 1.0/2*density*A*(v[uRow*nColsV + uCol] + v[uRow*nColsV + uCol - 1]);
			}
			

			// Explicit V-coefficients

			if (vRow == 0){
				// top
				fwV = 1.0/2*density*A*(u[vRow*nColsU + vCol] + u[(nRowsU - 1)*nColsU + vCol]);
				feV = 1.0/2*density*A*(u[vRow*nColsU + vCol + 1] + u[(nRowsU - 1)*nColsU + vCol + 1]);

				fsV = 1.0/2*density*A*(v[vRow*nColsV + vCol] + v[(vRow + 1)*nColsV + vCol]);
				fnV = 1.0/2*density*A*(v[vRow*nColsV + vCol] + v[(nRowsV - 1)*nColsV + vCol]);
			} else if(vRow == nRowsV - 1){
				// bottom
				fwV = 1.0/2*density*A*(u[(vRow - 1)*nColsU + vCol] + u[(0)*nColsU + vCol]);
				feV = 1.0/2*density*A*(u[(vRow - 1)*nColsU + vCol + 1] + u[(0)*nColsU + vCol + 1]);

				fsV = 1.0/2*density*A*(v[vRow*nColsV + vCol] + v[vCol]);
				fnV = 1.0/2*density*A*(v[vRow*nColsV + vCol] + v[(vRow - 1)*nColsV + vCol]);
			} else{
				// not a boundary
				fwV = 1.0/2*density*A*(u[vRow*nColsU + vCol] + u[(vRow - 1)*nColsU + vCol]);
				feV = 1.0/2*density*A*(u[vRow*nColsU + vCol + 1] + u[(vRow - 1)*nColsU + vCol + 1]);

				fsV = 1.0/2*density*A*(v[vRow*nColsV + vCol] + v[(vRow + 1)*nColsV + vCol]);
				fnV = 1.0/2*density*A*(v[vRow*nColsV + vCol] + v[(vRow - 1)*nColsV + vCol]);
			}

			// Check if U is in a solid interface. If not gather coefficients and calculate explicit component

			bool discFlagX = true;		// flag is created to avoid illegal and redundant memory accesses

			if(uCol == 0){
				if(Grid[uRow*nColsV + uCol] == 1){
					uExp[uRow*nColsU + uCol] = 0;
					uCoeff[uRow*nColsU + uCol] = 1;
					discFlagX = false;
				}
			} else if(uCol == nColsU - 1){
				if(Grid[uRow*nColsV + uCol - 1] == 1){
					uExp[uRow*nColsU + uCol] = 0;
					uCoeff[uRow*nColsU + uCol] = 1;
					discFlagX = false;
				}
			} else if(Grid[uRow*nColsV + uCol] == 1 || Grid[uRow*nColsV + uCol - 1] == 1){
				uExp[uRow*nColsU + uCol] = 0;
				uCoeff[uRow*nColsU + uCol] = 1;
				discFlagX = false;
			} 
			if(discFlagX == true){
				// Hybrid Discretization
				awU = findMax(fwU, d + fwU/2, 0.0f);
				aeU = findMax(-feU, d - feU/2, 0.0f);
				asU = findMax(fsU, d + fsU/2, 0.0f);
				anU = findMax(-fnU, d - fnU/2, 0.0f);

				DeltaF = feU - fwU + fnU - fsU;

				// printf("DeltaF = %f\n", DeltaF);

				// Store central coefficient

				// Should we select these coefficient better according to boundaries?

				uCoeff[uRow*nColsU + uCol] = awU + aeU + asU + anU + DeltaF;

				// Initalize uExplicit as a component of the previous iterative step

				uExp[uRow*nColsU + uCol] = (1.0f - o->alphaRelax)*u[uRow*nColsU + uCol];

				// Increment uExp based on neighborhood

				if(uCol != 0){
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(awU*u[uRow*nColsU + uCol - 1]);
				}

				if(uCol != nColsU - 1){
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(aeU*u[uRow*nColsU + uCol + 1]);
				}

				if(uRow != 0){
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(anU*u[(uRow - 1)*nColsU + uCol]);
				} else{
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(anU*u[(nRowsU - 1)*nColsU + uCol]);
				}

				if(uRow != nRowsU - 1){
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(asU*u[(uRow + 1)*nColsU + uCol]);
				} else{
					uExp[uRow*nColsU + uCol] += o->alphaRelax/uCoeff[uRow*nColsU + uCol]*(asU*u[(0)*nColsU + uCol]);
				}
			}

			// Check if V is in a solid interface. If not gather coefficients and calculate explicit component

			bool discFlag = true;			// flag is created to avoid illegal and redundant memory accesses

			if(vRow == 0){
				if(Grid[vRow*nColsV + vCol] == 1){
					// North boundary and grid is solid
					vExp[vRow*nColsV + vCol] = 0;
					vCoeff[vRow*nColsV + vCol] = 1;
					discFlag = false;
				}	
			}else if(vRow == nRowsV - 1){
				if(Grid[(vRow - 1)*nColsV + vCol] == 1){
					// South boundary and grid is solid
					vExp[vRow*nColsV + vCol] = 0;
					vCoeff[vRow*nColsV + vCol] = 1;
					discFlag = false;
				}
			} else if(Grid[vRow*nColsV + vCol] == 1 || Grid[(vRow - 1)*nColsV + vCol] == 1){
				vExp[vRow*nColsV + vCol] = 0;
				vCoeff[vRow*nColsV + vCol] = 1;
				discFlag = false;
			}

			if(discFlag == true){
				// Hybrid Discretization
				awV = findMax(fwV, d + fwV/2, 0.0f);
				aeV = findMax(-feV, d - feV/2, 0.0f);
				asV = findMax(fsV, d + fsV/2, 0.0f);
				anV = findMax(-fnV, d - fnV/2, 0.0f);

				DeltaF = feV - fwV + fnV - fsV;

				// Store central coefficient

				vCoeff[vRow*nColsV + vCol] = awV + aeV + asV + anV + DeltaF;

				// Initalize vExplicit as a component of the previous iterative step

				vExp[vRow*nColsV + vCol] = (1.0f - o->alphaRelax)*v[vRow*nColsV + vCol];

				// Increment vExp based on neighborhood

				if(vCol != 0){
					vExp[vRow*nColsV + vCol] += o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(awV*v[vRow*nColsV + vCol - 1]);
				}

				if(vCol != nColsV-1){
					vExp[vRow*nColsV + vCol] += o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(aeV*v[vRow*nColsV + vCol + 1]);
				}

				if(vRow != 0){
					vExp[vRow*nColsV + vCol] += o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(anV*v[(vRow - 1)*nColsV + vCol]);
				} else{
					vExp[vRow*nColsV + vCol] += o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(anV*v[(nRowsV - 1)*nColsV + vCol]);
				}

				if(vRow != nRowsV - 1){
					vExp[vRow*nColsV + vCol] += o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(asV*v[(vRow + 1)*nColsV + vCol]);
				} else{
					vExp[vRow*nColsV + vCol] += o->alphaRelax/vCoeff[vRow*nColsV + vCol]*(asV*v[(0)*nColsV + vCol]);
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
	memset(A, 0, sizeof(float)*info->numCellsX*info->numCellsY*5);

	// Indexing variables

	int nColsU = info->numCellsX + 1;
	int nColsV = info->numCellsX;
	int nColsP = info->numCellsX;
	int index = 0;

	// Neighborhood parameters

	float aE, aW, aS, aN, aP;
	float bp;

	for(int row = 0; row<info->numCellsY; row++){
		for(int col = 0; col<info->numCellsX; col++){
			index = row*nColsP + col;

			// Initialize to 0's to avoid NaN's
			aP = 0.0;
			aN = 0.0;
			aS = 0.0;
			aW = 0.0;
			aE = 0.0;
			bp = 0.0;
			RHS[index] = 0.0;

			// Check if P[index] is solid. If so, equation is P = 0 (central coeff = 1, RHS = 0)
			if(Grid[index] == 1){
				A[index*5 + 0] = -1.0;
				A[index*5 + 1] = 0.0;
				A[index*5 + 2] = 0.0;
				A[index*5 + 3] = 0.0;
				A[index*5 + 4] = 0.0;
				RHS[index] = 0.0;
			} else if(Grid[index] == 2)
			{
				A[index*5 + 0] = -1.0;
				A[index*5 + 1] = 0.0;
				A[index*5 + 2] = 0.0;
				A[index*5 + 3] = 0.0;
				A[index*5 + 4] = 0.0;
				RHS[index] = -o->PR;
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
					// aP += aW;
					aW = 0;
					// Check that east grid is not solid. If it is, set coeff to zero
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
					// aP += aE;
					aE = 0;
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

				if(col == 0){
					A[index*5 + 0] = -1.0;
					A[index*5 + 1] = 0.0;
					A[index*5 + 2] = 0.0;
					A[index*5 + 3] = 0.0;
					A[index*5 + 4] = 0.0;
					RHS[index] = -o->PL;
				} else if(col == nColsP - 1){
					A[index*5 + 0] = -1.0;
					A[index*5 + 1] = 0.0;
					A[index*5 + 2] = 0.0;
					A[index*5 + 3] = 0.0;
					A[index*5 + 4] = 0.0;
					RHS[index] = -o->PR;
				}
			}
		}
	}

	// Solve

	// Declare GPU arrays

	float *d_x_vec = NULL;
	float *d_temp_x_vec = NULL;
	
	float *d_Coeff = NULL;
	float *d_RHS = NULL;

	// Initialize the GPU arrays

	if(!initializeGPU(&d_x_vec, &d_temp_x_vec, &d_RHS, &d_Coeff, info))
	{
		printf("\n Error when allocating space in GPU");
		unInitializeGPU(&d_x_vec, &d_temp_x_vec, &d_RHS, &d_Coeff);
		return 0;
	}

	// Call the solver

	JacobiGPU2D(A, RHS, Pressure, o, info, d_x_vec, d_temp_x_vec, d_Coeff, d_RHS);

	// Un-initialize all dynamically allocated arrays

	unInitializeGPU(&d_x_vec, &d_temp_x_vec, &d_RHS, &d_Coeff);

	free(A);
	free(RHS);

	return 0;

}


int momentumCorrection(unsigned int *Grid, float *uExp, float *vExp, float* u, float* v,
	float *uCoeff, float *vCoeff, float *Pressure, options* o, simulationInfo* info)
{

	/*
	Momentum Correction:

	Inputs:
		- unsigned int *Grid: pointer to array holding boundary information
		- float *uExp: pointer to array to store the explicit component of u-velocity
		- float *vExp: pointer to array to store the explicit component of v-velocity
		- float *u: poiner to array storing the u-velocity
		- float *v: pointer to array storing the v-velocity
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
	float PL = o->PL;
	float PR = o->PR;

	// Useful indexing variables

	int nColsU = info->numCellsX + 1;
	int nColsV = info->numCellsY;
	int nColsP = info->numCellsY;

	int nRowsV = info->numCellsY+1;
	int nRowsP = info->numCellsY;

	int uRow, uCol, vRow, vCol;

	// Explicit Velocity Correction

	for(int i = 0; i<info->numCellsX+1; i++){
		for(int j = 0; j<info->numCellsY; j++){
			uRow = j;
			uCol = i;

			vRow = i;
			vCol = j;

			// Update U

			if(uCol == 0)
			{
				if(Grid[uRow*nColsP + uCol] == 1 || Grid[uRow*nColsP + uCol] == 2)
				{
					u[uRow*nColsU + uCol] = 0.0;
				} else
				{
					u[uRow*nColsU + uCol] = uExp[uRow*nColsU + uCol] + alpha*Area/uCoeff[uRow*nColsU + uCol]*(PL - Pressure[uRow*nColsP + uCol]);
				}
			}else if(uCol == info->numCellsX)
			{
				if(Grid[uRow*nColsP + uCol - 1] == 1 || Grid[uRow*nColsP + uCol - 1] == 2)
				{
					u[uRow*nColsU + uCol] = 0.0;
				}else
				{
					u[uRow*nColsU + uCol] = uExp[uRow*nColsU + uCol] + alpha*Area/uCoeff[uRow*nColsU + uCol]*(Pressure[uRow*nColsP + uCol - 1] - PR);
				}
			} else
			{
				if(Grid[uRow*nColsP + uCol] == 1 || Grid[uRow*nColsP + uCol - 1] == 1 || Grid[uRow*nColsP + uCol] == 2 || Grid[uRow*nColsP + uCol - 1] == 2)
				{
					u[uRow*nColsU + uCol] = 0.0;
				}else
				{
					u[uRow*nColsU + uCol] = uExp[uRow*nColsU + uCol] + alpha*Area/uCoeff[uRow*nColsU + uCol]*(Pressure[uRow*nColsP + uCol - 1] - Pressure[uRow*nColsP + uCol]);
				}
			}

			// Update V

			if(vRow == 0)
			{
				if(Grid[vRow*nColsV + vCol] == 1 || Grid[(nRowsP - 1)*nColsV + vCol] == 1 || Grid[vRow*nColsV + vCol] == 2 || Grid[(nRowsP - 1)*nColsV + vCol] == 2)
				{
					v[vRow*nColsV + vCol] = 0.0;
				} else
				{
					v[vRow*nColsV + vCol] = vExp[vRow*nColsV + vCol] + alpha*Area/vCoeff[vRow*nColsV + vCol]*(Pressure[vRow*nColsP + vCol] - Pressure[(nRowsP - 1)*nColsP + vCol]);
				}
			} else if(vRow == nRowsV - 1)
			{
				if(Grid[(vRow - 1)*nColsV + vCol] == 1 || Grid[vCol] == 1 || Grid[(vRow - 1)*nColsV + vCol] == 2 || Grid[vCol] == 2)
				{
					v[vRow*nColsV + vCol] = 0.0;
				} else
				{
					v[vRow*nColsV + vCol] = vExp[vRow*nColsV + vCol] + alpha*Area/vCoeff[vRow*nColsV + vCol]*(Pressure[vCol] - Pressure[(vRow - 1)*nColsP + vCol]);
				}
			} else
			{
				if(Grid[(vRow - 1)*nColsV + vCol] == 1 || Grid[(vRow)*nColsV + vCol] == 1 || Grid[(vRow - 1)*nColsV + vCol] == 2 || Grid[(vRow)*nColsV + vCol] == 2)
				{
					v[vRow*nColsV + vCol] = 0.0;
				}else
				{
					v[vRow*nColsV + vCol] = vExp[vRow*nColsV + vCol] + alpha*Area/vCoeff[vRow*nColsV + vCol]*(Pressure[vRow*nColsP + vCol] - Pressure[(vRow - 1)*nColsP + vCol]);
				}
			}
		}
	}

	return 0;

}