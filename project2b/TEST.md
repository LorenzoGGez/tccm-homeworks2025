In this file the main.f90 program has been run 4 different times. The aim
was comparing the results obtained by 4 different input file, with 4 
different matrix dimension. Here follows the table displating the results
af the calculations. In the first column the matrix dimension is insert, in the second
the time for CPU with fortran function (ff), in the third column the GPU time and
in the last one the CPU explicit function time (ef). Time data are in seconds.

n   CPU (ff)    GPU CPU (ef)
1000    0,7285  0,6783  7,6258E-03
2000    3,8167  1,1583  29,3038
4000    12,4187 0,9678  119,8204
8000    71,4385 2,4844  735,4859

