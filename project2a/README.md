# Sparse matrix multiplication program 

## ℹ️ About The Project

This repository contains a simple program that calculates the product of two simmetric matrices give in input in the sparse format (CSR). 
The product is obtained in three different way, the first subroutine exploits the property of the sparse matrices, using them in this format in order to obtain
directly the product. The second subroutine use an hand written code for the matrix product, using the rigorous mathematics formula,
after reconstructing the matrices in the non sparse-format. The third kind of calculation use the DGEMM routine from the BLAS library, that like the previous,
calculates the matrices product starting from the non-sparse format matrices. For each of this operation the cpu_time is obtained in order to understand the scaling
behaviour and how the filling rate influence the time of the calculations.

### Built With

![Fortran 90](https://img.shields.io/badge/Fortran_90-734F96?style=for-the-badge&logo=fortran&logoColor=white)

## 🚀 Getting Started

To run a simulation, you need to compile the `.f90` files using your preferred Fortran 90 compiler.
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

The program requires two input files (2 matrices), in the sparse format (CSR).

After launching the program, it will prompt you for the inputs filename, the print of the matrices and the number of iteration for the subroutine calls.
It will then display the the time that each routine needed in order to obtaine the final product.

Some test run on some function of the program can be achieved by compiling the test program.
Here we have an example on how to do it with gfortran:

```sh
    gfortran src/sparse_mod.f90 test/test.f90 -o testprog -lblas
```
