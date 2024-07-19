# PixelBasedPermeability
This repository is dedicated to the simulation of permeability in 2D structures via the Finite Volume Method (FVM). This approach was designed for maximum efficiency when generating large datasets for machine learning applications, thus uses the pixel resolution of the image as the base mesh for the simulation. Below is basic information on how to compile and run this code. For more detailed information about the code itself, refer to the documentation pdf. For more information regarding the computational model, refer to the publication (in preparation).

This repository includes three versions of the same code: a CPU, GPU, and HPC versions. The CPU version has an OpenMP-based accelerator for the iterative solver, and the user has the option of selecting the number of CPU cores. For the GPU version, the only user option is to use a CUDA capable GPU for acceleration. In this case, additional GPUs won't do anything. The HPC version can be used with multiple GPUs, and it requires at least one CPU core per GPU.

## Requirements

This list reflects what we tested on and can confirm that runs properly, but older versions might work.
- a
