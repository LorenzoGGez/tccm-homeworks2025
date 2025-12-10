 INSTALL instruction.

With this file we want to furnish the user with the minimal instruction for compiling and 
running the prrogram. Once installed the software, the user will be furnished of
the moldyn.f90 code. This is the code used for molecular dynamic calculation.
This code it's useless without an input file which contain all the data used for the
calculation. As a matter of fact, the code moldyn.f90, it's designed for
reading an external file and get from this external file all the necessary data for the calculation.
After the installation the user will also get a test file, which is a file with randomic
data, where the user will be able to understand how to set his own input data file, with
the right positions of the different data requested by the code. 
Once prepared the input data file, as shown in test file, the user needs to compile
the moldyn.f90 program with the following command: $ gfortran moldyn.f90 -o executable
The user can replace the word "executable" in the previous command with any word he likes.
At this point the program has been compiled, and the user will be provided with an executable
file, named as he chose before. The last thing to do now, is to launch the calculation of the executable.
For this pourpose, it's easly required to launch the following command: $ ./executable
Also in this case the "executable" is reaplaced by the word chosen at the previous compiling step.
If you followed precisely each step of this INSTALL instruction, you will be able to compute
molecular dynamic problems and obtaining, hopefully, important outcomes for your task.
