 INSTALL instruction.

With this file we want to furnish the user with the minimal instruction for compiling and 
running the prrogram. Once installed the software, the user will be furnished of
the necessary code in the src directory. This is the code used for calculating matrix product 
in three different ways and analyzing time performance.
This code it's useless without an input file which contain all the data used for the
calculation. As a matter of fact, the code sparse.f90, it's designed for
reading external files and get from these external files all the necessary data for the calculation.
After the installation the user will also get a test directory, where the user will be able 
to understand how some of the subroutine work and if they actually do work.
The user needs to compile the sparse.f90 program, the sparse_mod.f90 module and the lblas library with the following command: 
$  gfortran src/sparse_mod.f90 src/sparse.f90 -o executable -lblas
The user can replace the word "executable" in the previous command with any word he likes.
At this point the program has been compiled, and the user will be provided with an executable
file, named as he chose before. The last thing to do now, is to launch the calculation of the executable.
For this pourpose, it's easly required to launch the following command: $ ./executable
Also in this case the "executable" is reaplaced by the word chosen at the previous compiling step.
If you followed precisely each step of this INSTALL instruction, you will be able to compute
matrix products and obtaining the time for the operation.
