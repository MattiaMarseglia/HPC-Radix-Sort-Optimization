# Radix Sort Optimization Project (MPI, CUDA, OpenMP)

## Project Overview
This project focuses on parallelizing and evaluating the performance of the "Radix Sort" algorithm using CUDA, MPI, and OpenMP. Each section provides solutions to the problem of parallelizing the algorithm, analyzes results with different configurations, and evaluates performance, speedup, and efficiency.

## Table of Contents
1. [CUDA Implementation](#cuda-implementation)
2. [MPI Implementation](#mpi-implementation)
3. [OpenMP Implementation](#openmp-implementation)
4. [General Setup Instructions](#general-setup-instructions)
5. [License](#license)

## CUDA Implementation

### Problem Description
Parallelize and evaluate the performance of the Radix Sort algorithm using CUDA. Analyze results with different memory allocations (Global, Texture, Shared).

### Experimental Setup
- **Hardware**: GPU with Compute Capability (CC)
- **Software**: CUDA Toolkit, NVCC Compiler

### Report Requirements
- Problem description
- Experimental setup
- Hardware used (GPU CC)
- Software used
- Compiler configurations, options, workload descriptions, environment variables, etc.
- Performance metrics (Speedup, Efficiency)
- Result analysis and motivations
- Different solutions and optimal configuration
- Test case descriptions
- API descriptions (if any)
- Reproducibility instructions
- Colab notebook (if used, include as attachment)

### Source Code
- Provide makefile, source code, and folders to test and reproduce results.

## MPI Implementation

### Problem Description
Parallelize and evaluate the performance of the Radix Sort algorithm using MPI. Input data is provided from a file.

### Experimental Setup
- **Hardware**: CPUs, RAM, Virtual machines (if any), Colab setup (if any)
- **Software**: Operating system version, swap file/partition size, library versions, compiler versions, etc.

### Report Requirements
- Problem description
- Experimental setup
- Hardware used (CPUs, RAM, VM configuration, Colab setup)
- Software used
- Compiler configurations, options, workload descriptions, environment variables, etc.
- Performance metrics (Speedup, Efficiency)
- Result analysis and motivations
- Different solutions and optimal configuration
- Test case descriptions
- API descriptions (if any)
- Reproducibility instructions
- Colab notebook (if used, include as attachment)

### Source Code
- Provide makefile, source code, and folders to test and reproduce results.

## OpenMP Implementation

### Problem Description
Parallelize and evaluate the performance of the Radix Sort algorithm using OpenMP. Input data is provided from a file.

### Experimental Setup
- **Hardware**: CPUs, RAM, Virtual machines (if any), Colab setup (if any)
- **Software**: Operating system version, swap file/partition size, library versions, compiler versions, etc.

### Report Requirements
- Problem description
- Experimental setup
- Hardware used (CPUs, RAM, VM configuration, Colab setup)
- Software used
- Compiler configurations, options, workload descriptions, environment variables, etc.
- Performance metrics (Speedup, Efficiency)
- Result analysis and motivations
- Different solutions and optimal configuration
- Test case descriptions
- API descriptions (if any)
- Reproducibility instructions
- Colab notebook (if used, include as attachment)

### Source Code
- Provide makefile, source code, and folders to test and reproduce results.

## General Setup Instructions

### Setting Up the Environment
1. **CUDA**:
   - Install the CUDA Toolkit.
   - Ensure the NVCC compiler is properly configured.
   - Run `make` to build the CUDA implementation.

2. **MPI**:
   - Install MPI libraries (e.g., OpenMPI).
   - Compile the MPI code with `mpicc`.
   - Run the program using `mpirun`.

3. **OpenMP**:
   - Ensure your compiler supports OpenMP (e.g., GCC).
   - Compile the code with `gcc` using the `-fopenmp` flag.
   - Run the executable.

### Reproducing Results
- Detailed instructions on setting up the environment and running tests are provided in the respective sections of the report.

### License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

---

For more information, please refer to the respective implementation sections and their detailed reports.
