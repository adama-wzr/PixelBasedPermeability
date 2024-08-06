#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <stdbool.h>
#include <fstream>
#include <cfloat>
#include <set>
#include <string>

typedef struct VoronoiGenData
{
    int width;
    int height;
    int nseeds;
    int channelWidth;
    char filename[100];
    bool path;
    float porosity;
    int DistCalc;   // 1 for Euc, 2 for Manhattan
}VoronoiGenData;


typedef struct BatchOptsVor
{
    int nImg;
    int maxSeeds;
    int minSeeds;
    int maxChanR;
    int minChanR;
    int width;
    int height;
    int DistCalc;   // 1 for Euc, 2 for Manhattan, 3 for random
    bool pathFilter;    // true for filtering, false for not filtering
    int offset;
}BatchOptsVor;


typedef struct QSGSData
{
    int width;
    int height;
    float seedProb;
    char filename[100];
    bool path;
    float targetSVF;
    float DSide;
    float DDiag;
    float porosity;
}QSGSData;


typedef struct BatchOptsQSGS
{
    int nImg;
    float maxCd;
    float minCd;
    float maxDiag;
    float minDiag;
    float minR;
    float maxR;
    float minPs;
    float maxPs;
    int width;
    int height;
    bool pathFilter;    // true for filtering, false for not filtering
    int offset;
}BatchOptsQSGS;


typedef struct SphereData
{
    int width;
    int height;
    int sphereR;
    char filename[100];
    bool path;
    float targetSVF;
    float porosity;
}SphereData;


typedef struct BatchRSPM
{
    int nImg;
    int maxR;
    int minR;
    float maxSVF;
    float minSVF;
    bool pathFilter;
    int width;
    int height;
    int offset;
}BatchRSPM;

// define pair for coords

typedef std::pair<int, int> coordPair;

float calcPorosity2D(int *myImg, int height, int width)
{
    /*
        calcPorosity2D:

        Inputs:
            - *myImg: pointer to the image whose porosity is being calculated.
            - int height: height of the domain in pixels
            - int width: width of the domain in pixels
        Output:
            - None.

        Function will modify the variable "porosity" inside of the data structure.

    */

    float sum = 0;
    float totalPixels = width*height;

    for(int i = 0; i < width*height; i++)
    {
        if (myImg[i] == 0)
        {
            sum++;
        }
    }

    float pore = sum/totalPixels;

    return pore;
}


float EucDist(int i1, int j1, int i2, int j2){
    /*
        EucDist Function:

        Inputs:
            - four integers, locations of 2 points by coordinates.
        Output:
            - function returns the Euclidean distance between the two points as a float.

    */

    float dist = sqrt((float)(i2 - i1)*(i2 - i1) + (float)(j2 - j1)*(j2 - j1));

    return dist;
}

float ManhattanDist(int i1, int j1, int i2, int j2){
    /*
        ManhattanDist Function:

        Inputs:
            - four integers, locations of 2 points by coordinates.
        Output:
            - function returns the Manhattan distance between the two points as a float.

    */
    float dist = fabs((float)(i2 - i1)) + fabs((float)(j2 - j1));

    return dist;
}

int detectSides(unsigned char *arr, int height, int width, int ii, int jj){
    /*
        detectSides Function:

        Inputs:
            - pointer to arr: char array including morphological information from the images.
            - int height: height of the domain in pixels.
            - int width: width of the domain in pixels.
            - int ii and jj: coordinates, row and column (respectively), of the current pixel.
        Output:
            - None.
        
        Function returns a 0 if one of the sides are the same as the current pixel, returns 1 otherwise.

    */

    int N;
    int S;
    int E;
    int W;

    int currentPhase;

    int index = ii*width + jj;

    if(ii + 1 == height){
        N = 0;
        S = ii - 1;
    } else if(ii-1 < 0){
        N = ii + 1;
        S = height - 1;
    } else{
        N = ii + 1;
        S = ii - 1;
    }

    if(jj + 1 == width){
        W = 0;
        E = jj - 1;
    } else if(jj-1 < 0){
        W = jj + 1;
        E = width - 1;
    } else{
        W = jj + 1;
        E = jj - 1;
    }

    currentPhase = arr[index];

    if(currentPhase != arr[N*width + jj] || currentPhase != arr[S*width + jj] 
        || currentPhase != arr[ii*width + E] || currentPhase != arr[ii*width + W]){
        return 0;
    } else{
        return 1;
    }
    printf("Error?\n");
    return 1;
}

int detectDiag(unsigned char *arr, int height, int width, int ii, int jj){
    /*
        detectDiag Function:

        Inputs:
            - pointer to arr: char array including morphological information from the images.
            - int height: height of the domain in pixels.
            - int width: width of the domain in pixels.
            - int ii and jj: coordinates, row and column (respectively), of the current pixel.
        Output:
            - None.
        
        Function returns a 0 if one of the diagonals are the same as the current pixel, 1 otherwise.

    */

    int N;
    int S;
    int E;
    int W;

    int currentPhase;

    int index = ii*width + jj;

    if(ii + 1 == height){
        N = 0;
        S = ii - 1;
    } else if(ii-1 < 0){
        N = ii + 1;
        S = height - 1;
    } else{
        N = ii + 1;
        S = ii - 1;
    }

    if(jj + 1 == width){
        W = 0;
        E = jj - 1;
    } else if(jj-1 < 0){
        W = jj + 1;
        E = width - 1;
    } else{
        W = jj + 1;
        E = jj - 1;
    }

    currentPhase = arr[index];

    if(currentPhase != arr[N*width + W] || currentPhase != arr[N*width + E] 
        || currentPhase != arr[S*width + E] || currentPhase != arr[S*width + W]){
        return 0;
    } else{
        return 1;
    }
    printf("Error?\n");
    return 1;
}

bool FloodFillPath(int *Grid, int height, int width){
    /*
        FloodFillPath algorithm:

        Inputs:
            - pointer to Grid: array containing the location of solids and fluids
            - int height: height of the domain in pixels
            - int width: width of the domain in pixels

        Outputs:
            - None

        Function will modify the Grid array. If fluid is connected between the two boundaries,
        function will return a true, otherwise false.
    */

    // Declare variables
    int nElements = height*width;

    int* Domain = (int *)malloc(sizeof(int)*nElements);
    int index;

    // Initialize solid and fluid voxels

    for(int i = 0; i<nElements; i++){
        if(Grid[i] == 1)
        {
            Domain[i] = 1;  // solid
        }else
        {
            Domain[i] = -1; // fluid
        }
    }

    // Search is done from left to right
    // Initialize all fluid in the left boundary

    std::set<coordPair> cList;

    for(int row = 0; row<height; row++){
        if(Domain[row*width] == -1){
            Domain[row*width] = 0;
            cList.insert(std::make_pair(row, 0));
        }
    }

    // begin search

    while(!cList.empty()){
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
            tempRow = height - 1;
        } else{
            tempRow = row - 1;
        }

        // Update list if necessary

        if(Domain[tempRow*width + tempCol] == -1){
            Domain[tempRow*width + tempCol] = 0;
            cList.insert(std::make_pair(tempRow, tempCol));
        }

        // South

        tempCol = col;

        // check periodic boundary

        if(row == height - 1){
            tempRow = 0;
        } else{
            tempRow = row + 1;
        }

        // Update list if necessary

        if(Domain[tempRow*width + tempCol] == -1){
            Domain[tempRow*width + tempCol] = 0;
            cList.insert(std::make_pair(tempRow, tempCol));
        }

        // West

        if(col != 0){
            tempCol = col - 1;
            tempRow = row;

            if(Domain[tempRow*width + tempCol] == -1){
                Domain[tempRow*width + tempCol] = 0;
                cList.insert(std::make_pair(tempRow, tempCol));
            }
        }

        // East

        if(col != width - 1){
            tempCol = col + 1;
            tempRow = row;

            if(Domain[tempRow*width + tempCol] == -1){
                Domain[tempRow*width + tempCol] = 0;
                cList.insert(std::make_pair(tempRow, tempCol));
            }
        }
    }

    // Look for connected pixels in the right boundary

    for(int row = 0; row<height;row++){
        int col = width - 1;
        index = row*width + col;
        if (Domain[index] == 0){
            // connected pixels found
            free(Domain);
            return true;
        }
    }

    // no connected pixels found

    free(Domain);
    return false;  
}

int QSGS(QSGSData *myQSGS){
    /*
        QSGS:

        Inputs:
            - *QSGS: pointer to struct containing the QSGS parameters.
        Outputs:
            - none
        
        Simple implementation of the QSGS algorithm.
    */

    int height = myQSGS->height;
    int width = myQSGS->width;
    float Dside = myQSGS->DSide;
    float Ddiag = myQSGS->DDiag;
    float Ps = myQSGS->targetSVF;
    float Cd = myQSGS->seedProb;

    // dynamically allocate matrix for storing the image
    float solidCount = 0;
    int iterCount = 0;

    unsigned char *img = (unsigned char *)malloc(sizeof(char)*height*width);
    int *NeighborFlag = (int *)malloc(sizeof(int)*height*width);

    // initialize the arrays in memory to all zeroes

    memset(img, 0, sizeof(char)*height*width);
    memset(NeighborFlag, 0, sizeof(int)*height*width);

    // Step 1: random seeds

    for(int i = 0; i<width*height; i++){
        if((float)rand()/(float)RAND_MAX <= Cd)
        {
            img[i] = 255;
            solidCount++;
        }
    }

    // Step 2: Particle growth

    /*
        We will use the following flag scheme:
        
        0 -> no growth (either already solid or no contact with a solid)
        1 -> side is solid (Dside applied)
        2 -> diag is solid (Ddiag applied)

    */

    int index;
    while(solidCount/(height*width) < Ps)
    {
        iterCount++;
        for(int i = 0; i<height; i++){
            for(int j = 0; j<width; j++){
                index = i*width + j;
                // Apply appropriate flags
                if(img[index] == 0)
                {
                    if(detectSides(img, width, height, i, j) == 0)
                    {
                        NeighborFlag[index] = 1;
                    }else if(detectDiag(img, width, height, i, j) == 0)
                    {
                        NeighborFlag[index] = 2;
                    }else
                    {
                        NeighborFlag[index] = 0;
                    }
                }else
                {
                    NeighborFlag[index] = 0;
                }
            }
        }

        // Random particle growth based on assigned flags

        for(int i = 0; i<height; i++){
            for(int j = 0; j<width; j++){
                index = i*width + j;
                float myNum = (float)(rand())/(float)(RAND_MAX);
                if (NeighborFlag[index] == 1 && myNum < Dside)
                {
                    img[index] = 255;
                    solidCount++;
                }else if(NeighborFlag[index] == 2 && myNum < Ddiag)
                {
                    img[index] = 255;
                    solidCount++;
                }
            }
        }
    }

    free(NeighborFlag);

    // check for path and porosity

    int *Grid = (int *)malloc(sizeof(int)*height*width);

    memset(Grid, 0, sizeof(int)*width*height);

    int counter1, counter2;
    counter1 = 0;
    counter2 = 0;

    for(int i = 0; i < height*width; i++){
        if(img[i] > 0){
            Grid[i] = 1;
            counter1++;
        }else{
            Grid[i] = 0;
            counter2++;
        }
    }

    // myQSGS->porosity = calcPorosity2D(Grid, height, width);
    myQSGS ->porosity = 1.0 - solidCount/(height*width);
    myQSGS->path = FloodFillPath(Grid, height, width);

    free(Grid);

    stbi_write_jpg(myQSGS->filename, width, height, 1, img, 100);


    free(img);
    return 0;
}

int JFA_Vor(VoronoiGenData *myVor)
{
    /*
        JFA_Vor:

        Inputs:
            - *myVor: pointer to struct containing the Voronoi diagram parameters.
        Outputs:
            - none
        
        Simple implementation of the Jump Flooding Algorithm for Voronoi diagram approximation.
        Original publication:
        http://www.comp.nus.edu.sg/~tants/jfa.html
    */
    // unpack data structure
    int height = myVor->height;
    int width = myVor->width;
    int nseeds = myVor->nseeds;
    int channelWidth = myVor->channelWidth;

    // declare and define image arrays array

    int *img = (int *)malloc(sizeof(int)*height*width);

    memset(img, -1, sizeof(int)*height*width);

    int *NewImg = (int *)malloc(sizeof(int)*height*width);

    memset(NewImg, -1, sizeof(int)*height*width);

    // declare array to hold the row/col of each seed:

    int *seedLoc = (int *)malloc(sizeof(int)*nseeds*2);

    // Generate random seeds

    for(int i = 0; i<nseeds; i++){
        int tempi = rand() % width;         // col
        int tempj = rand() % height;        // row

        seedLoc[i*2 + 0] = tempj;
        seedLoc[i*2 + 1] = tempi;

        int index = tempj*width + tempi;

        img[index] = i + 1;
        NewImg[index] = i + 1;
    }

    // Perform JFA algorithm

    for(int a = width/2; a >= 1; a = a/2){
        // scan all pixels
        for(int i = 0; i<height; i++){
            for(int j = 0; j<width; j++){
                // expand central pixel to the domain
                for(int ir = i - a; ir <= i + a; ir++){
                    for(int jr = j - a; jr <= j + a; jr++){
                        // account for boundaries
                        int tempIr, tempJr;
                        if (ir < 0){
                            tempIr = 0;
                        }else if(ir > height - 1){
                            tempIr = height - 1;
                        }else{
                            tempIr = ir;
                        }

                        if(jr < 0){
                            tempJr = 0;
                        }else if(jr > width - 1){
                            tempJr = width - 1;
                        }else{
                            tempJr = jr;
                        }
                            // img with index ir and jr never exists! fix all those indexes
                        if(img[i*width + j] == -1 && img[tempIr*width + tempJr] != -1){
                            img[i*width + j] = img[tempIr*width + tempJr];
                        }else if(img[i*width + j] != img[tempIr*width + tempJr] && img[tempIr*width + tempJr] != -1 && img[i*width + j] != -1){
                            float dist, distNB;
                            if(myVor->DistCalc == 1){
                                dist = EucDist(i, j, seedLoc[(img[i*width + j] - 1)*2 + 0], seedLoc[(img[i*width + j] - 1)*2 + 1]);
                                distNB = EucDist(i, j, seedLoc[(img[tempIr*width + tempJr] - 1)*2 + 0], seedLoc[(img[tempIr*width + tempJr] - 1)*2 + 1]);
                            }else if(myVor->DistCalc == 2){
                                dist = ManhattanDist(i, j, seedLoc[(img[i*width + j] - 1)*2 + 0], seedLoc[(img[i*width + j] - 1)*2 + 1]);
                                distNB = ManhattanDist(i, j, seedLoc[(img[tempIr*width + tempJr] - 1)*2 + 0], seedLoc[(img[tempIr*width + tempJr] - 1)*2 + 1]);
                            }
                            
                            if(dist > distNB){
                                img[i*width + j] = img[tempIr*width + tempJr];
                            }
                        }
                    }
                }
            }
        }
    }

    // Need to find pixels that are in the interface

    for(int i = 0; i<height; i++){
        for(int j = 0; j<width; j++){
            int a = 1;
            for(int ir = i - a; ir <= i + a; ir++){
                for(int jr = j - a; jr <= j + a; jr++){
                    // account for boundaries
                    int tempIr, tempJr;
                    if (ir < 0){
                        tempIr = 0;
                    }else if(ir > height - 1){
                        tempIr = height - 1;
                    }else{
                        tempIr = ir;
                    }

                    if(jr < 0){
                        tempJr = 0;
                    }else if(jr > width - 1){
                        tempJr = width - 1;
                    }else{
                        tempJr = jr;
                    }
                    if(img[i*width + j] != img[tempIr*width + tempJr] && img[i*width + j] != 0 && img[tempIr*width + tempJr] != 0)
                    {
                        img[tempIr*width + tempJr] = 0;
                    }
                }
            }
        }
    }

    // Expand channels

    memcpy(NewImg, img, sizeof(int)*width*height);       // dest, src, size

    int radius = channelWidth/2;

    for(int j = 0; j<height; j++){
        for(int i = 0; i<width; i++){
            if(img[j*width+i] == 0){
                for(int ir = i-radius; ir < i+radius+1; ir++){
                    for(int jr = j-radius; jr < j+radius+1; jr++){
                        int temp_i, temp_j;

                        if(ir < 0){
                            temp_i = 0;
                        } else if(ir >= width){
                            temp_i = width - 1;
                        } else{
                            temp_i = ir;
                        }

                        if(jr < 0){
                            temp_j = 0;
                        } else if(jr >= height){
                            temp_j = height - 1;
                        } else{
                            temp_j = jr;
                        }

                        if((ir - i)*(ir - i) + (jr - j)*(jr-j) <= radius*radius){
                            NewImg[temp_j*width + temp_i] = 0;
                        }
                    }
                }
            }
        }
    }

    // Copy and free memory

    memcpy(img, NewImg, sizeof(int)*width*height);       // dest, src, size

    free(NewImg);

    // Calculate porosity

    myVor->porosity = calcPorosity2D(img, height, width);

    // Run pathfinding algorithm

    int *pathGrid = (int *)malloc(sizeof(int)*width*height);

    for(int j = 0; j<height; j++){
        for(int i = 0; i<width; i++){
            if( img[j*width + i] > 0)
            {
                pathGrid[j*width + i] = 1;
            } else
            {
                pathGrid[j*width + i] = 0;
            }
        }
    }

    myVor->path = FloodFillPath(pathGrid, height, width);
    free(pathGrid);

    // Save image

    char *gen_img = (char*)malloc(sizeof(char)*width*height);

    for(int j = 0; j<height; j++){
        for(int i = 0; i<width; i++){
            if(img[j*width + i] > 0){
                gen_img[j*width + i] = 255;
            } else{
                gen_img[j*width + i] = 0;
            }
        }
    }
    
    stbi_write_jpg(myVor->filename, width, height, 1, gen_img, 100);

    free(gen_img);
    free(img);
   
   
   return 0;
}

int RSPM(SphereData *mySphere)
{
    /*
        Random Sphere Packing Method (RSPM)

        Inputs:
            - *mySphere = user entered data for RSPM

        Outputs:
            - none.

        Function will randomly pack a domain with overlapping spheres until the domain 
        reaches the target SVF. The image is then saved.
    */

    int height = mySphere->height;
    int width = mySphere->width;
    int radius = mySphere->sphereR;
    float SVF = 0.0;

    // declare and define array where img will be stored

    unsigned char* img = (unsigned char*)malloc(sizeof(char)*height*width);

    memset(img, 0, sizeof(char)*height*width);

    while(SVF < mySphere->targetSVF)
    {
        int centerRow = rand() % height;
        int centerCol = rand() % width;

        for(int i = centerRow - radius; i <= centerRow + radius; i++){
            for(int j = centerCol - radius; j <= centerCol + radius; j++){
                // Proceed only if within a circle
                if(EucDist(centerRow, centerCol, i, j) < radius){
                    int temp_i, temp_j;
                    if(i < 0){
                        temp_i = height + i;
                    }else if(i > height - 1){
                        temp_i = i - height;
                    }else{
                        temp_i = i;
                    }

                    if(j < 0){
                        temp_j = width + j;
                    } else if(j > width - 1){
                        temp_j = j - width;
                    }else{
                        temp_j = j;
                    }

                    img[temp_i*width + temp_j] = 255;
                }
            }
        }

        // Calculate SVF after drawing
        float solidCount = 0;
        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++){
                if(img[i*width + j] > 0){
                    solidCount++;
                }
            }
        }

        SVF = solidCount/(float)(height*width);

    }

    // Calculate porosity and find path

    int* Grid = (int *)malloc(sizeof(int)*width*height);

    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            if(img[i*width + j] > 0){
                Grid[i*width + j] = 1;
            }else{
                Grid[i*width + j] = 0;
            }
        }
    }

    mySphere->porosity = calcPorosity2D(Grid, height, width);
    mySphere->path = FloodFillPath(Grid, height, width);

    free(Grid);

    printf("Porosity = %f, Path = %d\n", mySphere->porosity, mySphere->path);

    // Save image

    stbi_write_jpg(mySphere->filename, width, height, 1, img, 100);

    free(img);

    return 0;
}

int BatchVoronoi(BatchOptsVor *opts)
{
    /*
        Batch Voronoi Tesselations:

        Inputs:
            - *opts = user entered options for batch generation of voronoi tessellations

        Outputs:
            - none.

        Function will generate the prescribed number of tessellations and save them.
    */

    // Declare and define useful arrays, structs

    VoronoiGenData *genData;

    genData = (VoronoiGenData *)malloc(sizeof(VoronoiGenData)*opts->nImg);

    // Begin Generation

    for(int n = 0; n < opts->nImg; n++)
    {
        // seed random number generator
        srand(time(NULL)*(n + 1));

        // set values randomly where applicable

        genData[n].width = opts->width;
        genData[n].height = opts->height;
        genData[n].nseeds = rand() % (opts->maxSeeds - opts->minSeeds) + opts->minSeeds;
        genData[n].channelWidth = rand() % (opts->maxChanR - opts->minChanR) + opts->minChanR;
        if(opts->DistCalc == 3){
            genData[n].DistCalc = rand() % 1 + 1;
        }else{
            genData[n].DistCalc = opts->DistCalc;
        }
        

        sprintf(genData[n].filename, "%05d.jpg", n+opts->offset);

        // Voronoi2D(&genData[n]);
        JFA_Vor(&genData[n]);
        int attemptCount = 1;
        if(opts->pathFilter == true && genData[n].path == 0)    // check path flag, and re-do image if necessary
        {
            while(genData[n].path == 0)
            {
                attemptCount++;
                genData[n].nseeds = rand() % (opts->maxSeeds - opts->minSeeds) + opts->minSeeds;
                genData[n].channelWidth = rand() % (opts->maxChanR - opts->minChanR) + opts->minChanR;
                if(opts->DistCalc == 3){
                    genData[n].DistCalc = rand() % 1 + 1;
                }else{
                    genData[n].DistCalc = opts->DistCalc;
                }
                JFA_Vor(&genData[n]);
                printf("Image num = %d, Attempt Count = %d\n", n, attemptCount);
            }
        }

    }

    // Create file with statistics

    FILE* VOR;

    VOR = fopen("VoronoiData.csv", "a+");

    fprintf(VOR, "img,pore,chanr,nseeds,path\n");

    for(int n = 0; n<opts->nImg; n++){
        fprintf(VOR, "%d,%f,%d,%d,%d\n", n, genData[n].porosity, genData[n].channelWidth, genData[n].nseeds, genData[n].path);
    }

    fclose(VOR);

    free(genData);

    return 0;
}

int BatchQSGS(BatchOptsQSGS *opts){
    /*
        Batch Quartet Structure Generation Set:

        Inputs:
            - *opts = user entered options for batch generation of QSGS

        Outputs:
            - none.

        Function will generate the prescribed number of QSGS and save them.
    */

   QSGSData  *genData;

   genData = (QSGSData *)malloc(sizeof(QSGSData)*opts->nImg);

   float scale;
   float R;

    for(int n = 0; n<opts->nImg; n++){
        // re-seed every iteration
        srand((n + 1)*time(NULL));

        // define and store essential variables for QSGS
        genData[n].width = opts->width;
        genData[n].height = opts->height;
        scale = rand()/(float) RAND_MAX;
        genData[n].seedProb = opts->minCd + scale*(opts->maxCd - opts->minCd);
        scale = rand()/(float) RAND_MAX;
        genData[n].targetSVF = opts->minPs + scale*(opts->maxPs - opts->minPs);
        scale = rand()/(float) RAND_MAX;
        genData[n].DDiag = opts->minDiag + scale*(opts->maxDiag - opts->minDiag);
        scale = rand()/(float) RAND_MAX;
        R = opts->minR + scale*(opts->maxR - opts->minR);
        genData[n].DSide = genData[n].DDiag*R;

        // define filename

        sprintf(genData[n].filename, "%05d.jpg", n+opts->offset);

        // run QSGS

        QSGS(&genData[n]);

        // check path flag and re-do image if necessary

        int attemptCount = 1;

        if(opts->pathFilter == true && genData[n].path == 0)    // check path flag, and re-do image if necessary
        {
            while(genData[n].path == 0)
            {
                attemptCount++;
                scale = rand()/(float) RAND_MAX;
                genData[n].seedProb = opts->minCd + scale*(opts->maxCd - opts->minCd);
                scale = rand()/(float) RAND_MAX;
                genData[n].targetSVF = opts->minPs + scale*(opts->maxPs - opts->minPs);
                scale = rand()/(float) RAND_MAX;
                genData[n].DDiag = opts->minDiag + scale*(opts->maxDiag - opts->minDiag);
                scale = rand()/(float) RAND_MAX;
                R = opts->minR + scale*(opts->maxR - opts->minR);
                genData[n].DSide = genData[n].DDiag*R;
                QSGS(&genData[n]);
                printf("Image num = %d, Attempt Count = %d\n", n, attemptCount);
            }
        }

    }

    // Create file with statistics

    FILE* QSGS_FILE;

    QSGS_FILE = fopen("QSGS_Data.csv", "a+");

    fprintf(QSGS_FILE, "img,pore,Cd,Ddiag,Dside,path\n");

    for(int n = 0; n<opts->nImg; n++){
        fprintf(QSGS_FILE, "%d,%f,%f,%f,%f,%d\n", n, genData[n].porosity, genData[n].seedProb, genData[n].DDiag, genData[n].DSide, genData[n].path);
    }

    fclose(QSGS_FILE);

    free(genData);

    return 0;
}

int BatchSphere(BatchRSPM *opts)
{
    /*
        Function BatchSphere

        Inputs:
            - pointer to BatchRSPM: pointer to BatchRSPM struct cointaining the batch options.
        Outputs: 
            - None.
        
        Function will generate random structures randomly packed with circles, based on user entered data.
    */

    // Declare and define datastructure to hold the data for n images

    SphereData *genData;
    genData = (SphereData *)malloc(sizeof(SphereData)*opts->nImg);
    
    // Declare useful variables
    float scale;

    for(int n = 0; n < opts->nImg; n++){
        // re-seed every iteration
        srand((n + 1)*time(NULL));

        genData[n].width = opts->width;
        genData[n].height = opts->height;
        scale = rand()/(float) RAND_MAX;
        genData[n].targetSVF = opts->minSVF + scale*(opts->maxSVF - opts->minSVF);
        genData[n].sphereR = rand() % (opts->maxR - opts->minR) + opts->minR;

        printf("Here\n");

        printf("r = %d, SVF = %f\n", genData[n].sphereR, genData[n].targetSVF);

        // define filename

        sprintf(genData[n].filename, "%05d.jpg", n+opts->offset);

        // Run Sphere Packing

        RSPM(&genData[n]);

        int attemptCount = 1;

        // Generate again if path is not found

        while(genData[n].path == 0){
            scale = rand()/(float) RAND_MAX;
            genData[n].targetSVF = opts->minSVF + scale*(opts->maxSVF - opts->minSVF);
            genData[n].sphereR = rand() % (opts->maxR - opts->minR) + opts->minR;

            RSPM(&genData[n]);
            printf("Image num = %d, Attempt Count = %d\n", n, attemptCount);
        }
    }

    // Create file with statistics

    FILE* RSPM_FILE;

    RSPM_FILE = fopen("Sphere_Data.csv", "a+");

    fprintf(RSPM_FILE, "img,pore,R,path\n");

    for(int n = 0; n<opts->nImg; n++){
        fprintf(RSPM_FILE, "%d,%f,%d,%d\n", n, genData[n].porosity, genData[n].sphereR, genData[n].path);
    }

    fclose(RSPM_FILE);

    free(genData);

    return 0;
}