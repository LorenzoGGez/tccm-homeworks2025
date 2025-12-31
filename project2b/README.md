The program purchased is usefull to compute the first eigenvalue of a symmetric matrix.
This purpose is accomplished by two CPU implementation, using fortran function
and using the non parallelized GPU code, and with a GPU implementation using OpenMP.
The program needs an input file. This file shall be named 'input' and have in it's first line
the matrix dimension, just a single number since the matrix is symmetry. In the second line 
the convergence threshold shall be written while on third line you shall write the maximum
iterations number. 
