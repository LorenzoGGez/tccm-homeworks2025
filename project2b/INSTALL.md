Once purchased this program it's usage it's quite simple.
Use the following command to compile the program in the directory
where it is downloaded:
nvfortran -mp=gpu -gpu=cc70 main.f90 -o main 
Keep in mind that this command works once loaded nvhpc module.
In order to run the calculation, once built an input file (see README.md),
the calculation starts by submitting sbatch.sh script. If you followed
precisely these instructions, you would be fine with your calculation, which
will display an output file with your results. Do not forget that this 
program works on CINECA super computer.
