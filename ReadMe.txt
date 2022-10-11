C language algorithm for matrix-matrix product calculus in a parallel environment
for distributed memory MIMD architectures using MPI library and BMR communication strategy.
Two matrixes are distributed to pxp processes, which are arranged according to a
two-dimensional grid topology. 
This algorithm implements a BMR communication strategy and it's organized to use a subroutine
for the construction of the two-dimensional process grid and a subroutine for the calculation
of local products.
