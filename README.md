# PixelBasedPermeability
This repository is dedicated to the simulation of permeability in 2D structures via the Finite Volume Method (FVM). This approach was designed for maximum efficiency when generating large datasets for machine learning applications, thus uses the pixel resolution of the image as the base mesh for the simulation. Below is basic information on how to compile and run this code. For more detailed information about the code itself, refer to the documentation pdf. For more information regarding the computational model, refer to the publication (in preparation).

This repository includes three versions of the same code: a CPU, GPU, and HPC versions. The CPU version has an OpenMP-based accelerator for the iterative solver, and the user has the option of selecting the number of CPU cores. For the GPU version, the only user option is to use a CUDA capable GPU for acceleration. In this case, additional GPUs won't do anything. The HPC version can be used with multiple GPUs, and it requires at least one CPU core per GPU.

# Table of Contents

1. [Requirements](#requirements)
2. [CPU Compilation](#cpu-compilation)
3. [GPU Compilation](#gpu-compilation)
4. [HPC Compilation](#hpc-compilation)
5. [Required Files](#required-files)
6. [Data Generation](#data-generation)
7. [How to Cite](#how-to-cite)
8. [Authors](#authors)
9. [Documentation](#documentation)
10. [Acknowledgements](#acknowledgements)

## Requirements

This list reflects what we tested on and can confirm that runs properly, but older versions might work. Might work with other compilers as well.
- NVIDIA Compute capability >= 8.6
- CUDA >= 11.0
- gcc >= 11.0
- OpenMP >= 4.5
- [stb_image](https://github.com/nothings/stb) latest version

The code has been tested on Ubuntu >= 20.04, Windows 10 and 11. The HPC version has only been tested on Rocky Linux 8.7.

## CPU Compilation

Need to invoke the OpenMP from gcc. On Ubuntu, assuming all files are in the same folder:

```bash
g++ -fopenmp Perm2D.cpp
```
## GPU Compilation

With the NVIDIA suite installed properly and already added to the path, also assuming all required files are in the same folder.

```bash
nvcc Perm2D.cu
```
## HPC Compilation

The HPC compilation will depend on the HPC environment. Assuming all modules are loaded properly and the appropriate files are in the same folder.

```bash
nvcc -Xcompiler -fopenmp Perm2D.cu
```
## Required Files

All these files have to be in the same folder (or in the path for compilation/run).

- 2D grayscale .jpg image.
- Main Perm2D file (.cpp or .cu)
- Helper Perm2D file (.h or .cuh)
- input.txt
- stb_image.h

## Data Generation

The folder [2D Image Generation](https://github.com/adama-wzr/PixelBasedPermeability/tree/main/2D%20Image%20Generation) is included as it is an important part of the publication. While that code can be used to reproduce the dataset of the publication, it is also a great source of data for any user that is looking to test the code or even generate a dataset for something else. The code is capable of generating batches of QSGS, Voronoi, and Random-Sphere Packing Method, all in 2D, with relative efficiency.

## How to Cite

Publication is in preparation at the moment. If you need to use this code and there is no publication available yet, contact one of the authors. There will be one publication for the code, and one link to the open source dataset, which was generated using this code.

## Authors

- Main developer: Andre Adam (The University of Kansas)
    - [ResearchGate](https://www.researchgate.net/profile/Andre-Adam-2)
    - [GoogleScholar](https://scholar.google.com/citations?hl=en&user=aP_rDkMAAAAJ)
    - [GitHub](https://github.com/adama-wzr)
- Advisor: Dr. Xianglin Li (Washingtion University in St. Louis)
    - [Website](https://xianglinli.wixsite.com/mysite)
    - [GoogleScholar](https://scholar.google.com/citations?user=8y0Vd8cAAAAJ&hl=en)
- Advisor: Dr. Huazhen Fang (The University of Kansas)
    - [Website](https://fang.ku.edu/)
    - [Lab Website](https://www.issl.space/)
    - [GoogleScholar](https://scholar.google.com/citations?user=3m7Yd4YAAAAJ&hl=en)
- Validation: Silven Stallard (The University of Kansas)
    - [ResearchGate](https://www.researchgate.net/profile/Silven_Stallard)

## Documentation

The publication (upcoming) is an excellent source of basic information on the formulation and validation. The documentation pdf is a more in-depth source on the mathematical formulation and code implementation, while also providing technical insight on how to run and modify the code included in this repository. The documentation for the image generation and permeability code are separate documents, but both of them can be found [here](https://github.com/adama-wzr/PixelBasedPermeability/tree/main/Documentation).

## Acknowledgements

This work wouldn't be possible without the computational time awarded as part of the following grants:

This work used Expanse(GPU)/Bridges2(CPU) at SDSC/PSC through allocations MAT210014 and MAT230071 from the Advanced Cyberinfrastructure Coordination Ecosystem: Services & Support (ACCESS) program, which is supported by National Science Foundation grants #2138259, #2138286, #2138307, #2137603, and #2138296.

X.L. highly appreciates the support from the National Science Foundation (Award 1941083 and 2329821), A.A. and S.S. appreciate the funding support by NASA EPSCoR (Award 80NSSC22M0221).

