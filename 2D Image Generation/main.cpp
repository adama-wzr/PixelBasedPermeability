#include "helper.h"

int main(void){

    bool batchFlag, pathFilter;

    batchFlag = true;
    pathFilter = true;

    /*
        For standalone 2D image generation, use the following:
        - 0 for Voronoi
        - 1 for QSGS
        - 2 for Circles
    */

    int genMethod = 0;

    // Offset for filename numbering

    int offset = 0;

    // Good practice to write data properly

    fflush(stdout);

    if(batchFlag == false)
    {

        // seed random number generator

        srand(time(NULL));
        if(genMethod == 0)
        {
            // Declare and define all parameters necessary
            VoronoiGenData myVor;
            myVor.channelWidth = 4;
            myVor.height = 128;
            myVor.width = 128;
            myVor.nseeds = 8;
            myVor.DistCalc = 1;
            sprintf(myVor.filename, "VoronoiTest.jpg");
            // Generate Image
            JFA_Vor(&myVor);
        }else if(genMethod == 1)
        {
            // Declare and define all parameters necessary
            QSGSData myQSGS;
            myQSGS.width = 128;
            myQSGS.height = 128;
            myQSGS.seedProb= 0.0025;
            sprintf(myQSGS.filename, "QSGSTest.jpg");
            myQSGS.targetSVF = 0.75;
            myQSGS.DSide = 0.0004;
            myQSGS.DDiag = 0.0001;
            // Generate Image
            QSGS(&myQSGS);
        }else if(genMethod == 2)
        {
            // Declare and define all parameters necessary
            SphereData mySphere;
            mySphere.width = 128;
            mySphere.height = 128;
            mySphere.targetSVF = 0.8;
            mySphere.sphereR = 10;
            sprintf(mySphere.filename, "SphereTest.jpg");
            // Generate image
            RSPM(&mySphere);
        }else
        {
            printf("Error: Method entered not available.\n");
            return 1;
        }  
    } else{
        if(genMethod == 0)
        {
            // Declare batch options
            BatchOptsVor opts;

            // Define options

            opts.nImg = 300;
            opts.maxSeeds = 64;
            opts.minSeeds = 4;
            opts.maxChanR = 12;
            opts.minChanR = 2;
            opts.width = 128;
            opts.height = 128;
            opts.DistCalc = 3;
            opts.pathFilter = pathFilter;
            opts.offset = offset;

            BatchVoronoi(&opts);
        }else if(genMethod == 1)
        {
            // Declare batch options
            BatchOptsQSGS opts;
            
            // Define options
            opts.nImg = 300;
            opts.maxCd = 0.005;
            opts.minCd = 0.0001;
            opts.maxDiag = 0.0004;
            opts.minDiag = 0.0001;
            opts.minR = 1;
            opts.maxR = 16;
            opts.minPs = 0.03;
            opts.maxPs = 0.8;
            opts.width = 128;
            opts.height = 128;
            opts.pathFilter = pathFilter;
            opts.offset = offset;

            // Call batch function

            BatchQSGS(&opts);
        }else if(genMethod == 2)
        {
            // Declare Batch options
            BatchRSPM opts;
            // Define options
            opts.nImg = 300;
            opts.pathFilter = pathFilter;
            opts.height = 128;
            opts.width = 128;
            opts.maxR = 20;
            opts.minR = 2;
            opts.maxSVF = 0.9;
            opts.minSVF = 0.05;
            opts.offset = offset;

            // Call Batch
            BatchSphere(&opts);

        }else
        {
            printf("Error. Method entered not available.\n");
            return 1;
        }
        
    }

    


    return 0;
}