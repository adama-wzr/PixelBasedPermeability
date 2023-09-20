#include "Perm2D.cuh"

int main(void){

	// More efficient printing with parallel computing/Linux
	fflush(stdout);

	// Parse user entered options

	options opts;	// struct to hold options

	char inputFilename[100];

	sprintf(inputFilename, "input.txt");

	readInputFile(inputFilename, &opts);



	return 0;
}