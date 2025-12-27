# Noble Gas Cluster Molecular Dynamics Program

## ℹ️ About The Project

This repository contains a simple program that calculates the dynamics of a system of N atoms interacting under the Lennard-Jones potential.
The project was written as an exercise to learn how to build scientific software using Fortran 90.
We provide the `.f90` source code so that you can compile and use it on any operating system.

### Built With

![Fortran 90](https://img.shields.io/badge/Fortran_90-734F96?style=for-the-badge&logo=fortran&logoColor=white)

## 🚀 Getting Started

To run a simulation, you need to compile the `.f90` file using your preferred Fortran 90 compiler.
Below is a list of the most common ones:

| Compiler | Type | Supported OS | Best For |
| :--- | :--- | :--- | :--- |
| **[GNU Fortran (gfortran)](https://gcc.gnu.org/wiki/GFortran)** | 🟢 Open Source | 🐧 Linux<br>🪟 Windows<br>🍎 macOS | **General Use.** The standard for open source projects. |
| **[Intel OneAPI (ifx / ifort)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html)** | 🏢 Free (Proprietary) | 🐧 Linux<br>🪟 Windows<br>🍎 macOS | **Performance.** Highly optimized for Intel CPUs. |
| **[NVIDIA HPC SDK (nvfortran)](https://developer.nvidia.com/hpc-sdk)** | 🏢 Free (Proprietary) | 🐧 Linux<br>🪟 Windows (WSL) | **GPU Acceleration.** Best for CUDA/OpenACC. |
| **[LFortran](https://lfortran.org/)** | 🟢 Open Source | 🐧 Linux<br>🪟 Windows<br>🍎 macOS | **Interactive.** Modern, fast, and interactive usage. |

Once compiled, you will obtain an executable file to run the program.

### Installation

1.  Clone the repo
    ```sh
    git clone https://github.com/LorenzoGGez/tccm-homeworks2025.git
    ```

## 💻 Usage

The program requires an input file structured exactly like **[inp.txt](inp.txt)**; detailed instructions are provided inside the file itself.

After launching the program, it will prompt you for the input filename and the number of simulation steps.
It will then display the initial conditions and the final energy state on the screen, while generating a trajectory output file in the background.

Visualization is supported via **[Molden](https://www3.cmbi.umcn.nl/molden/)**.
