## PERFORMANCE ANALYSIS
In this file we are gonna discuss the behaviour of the three routines for matrices product in relation to dimension and the sparsity fill percentage of the matrices.
All the calculation were submitted to the same local machine and the times are taken for 100 calling of each subroutine.
First we are gonna look how the times are related to matrix dimension with a fix percentage of filling, then we will have a look to the dimension of the matrix 
and how. 


# MATRIX DIMENSION
With a fixed filling percentage we can compare the performance of the three routine in relation with the dimension increasing. 
For each of the algorithm we can clearly see that the increasing dimension lead to higher CPU time more or less in the same way for each ones. 
This mean that the dimension act on the three routine in the same way.
| Dimension (N) | Sparse Time (s) | Fortran Time (s) | DGEMM Time (s) |
|--------------:|----------------:|-----------------:|---------------:|
|            12 |        0.000059 |         0.000435 |       0.000123 |
|            25 |        0.000330 |         0.003990 |       0.001590 |
|            50 |        0.003911 |         0.042892 |       0.010769 |
|            75 |        0.011500 |         0.236512 |       0.053399 |
|           125 |        0.054276 |         1.216069 |       0.226622 |
|           250 |        0.571784 |         7.611464 |       1.868338 |
|           500 |        6.148291 |        78.596401 |      16.567160 |
|           750 |       20.206409 |       270.617534 |      50.661860 |
|           999 |       48.232266 |       448.086570 |      73.424045 |


<img width="1000" height="600" alt="spee vs dimensio" src="https://github.com/user-attachments/assets/2c0faf1c-145e-401a-9645-08eb9aaa6aaa" />



# MATRIX FILLING PERCENTAGE
Now we have a look on the performance behaviour with fixed dimension (999 X 999). Just the plot here is enough to understand the behaviour, 
the hand-written subroutine and the DGEMM routine are indipendent from the filling percentage of the matrices whis for the sparse algorithm 
we can see how dramtic are the results. In fact this algorithm is based on the sparse format of the matrices and we can see that for really sparse matrices 
the result are better than both the other algorithm. If we compare with the hand-written one we can see that in the 5%-10% range there 
is the switch in which one works better. With the DGEMM subroutine this range gets even lower (2%-5%), proving once more that the sparse algorithm it's good only for 
really sparse matrices.
| Density (%) | Op. Count (N_OPS) | Sparse Time (s) | Fortran Time (s) | DGEMM Time (s) |
| :---: | :---: | :---: | :---: | :---: |
| **1%** | 61,346 | 13.6120 | 447.8924 | 73.0426 |
| **2%** | 220,087 | 48.2323 | 448.0866 | 73.4240 |
| **5%** | 1,288,351 | 280.3684 | 427.7165 | 73.1279 |
| **10%** | 5,059,664 | 1,059.3003 | 616.7415 | 93.3578 |
| **25%** | 31,294,517 | 6,439.5221 | 429.0154 | 70.3953 |
| **50%** | 124,475,803 | 25,450.6417 | 702.4803 | 98.4092 |


<img width="1000" height="600" alt="speed_vs _density" src="https://github.com/user-attachments/assets/914a205c-25a5-4e28-a1dc-4d7907fe1b41" />


#NOTE
The previous data have some error due the usage of the machine with other operation, while doing these operation,
that may lead to increase the time for some of the subroutine.
