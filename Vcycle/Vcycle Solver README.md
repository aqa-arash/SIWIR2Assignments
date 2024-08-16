**V-Cycle Solver README**

**Introduction**

The program provides a geometric multigrid (GMG) solver that computes an approximate solution to a discretization of elliptic PDE using the V-cycle multigrid method. Additionally, this program utilizes OpenMP for parallelization and includes functionalities for direct solving, restriction, prolongation, and error calculation.

The elliptic partial differential equation is defined as,

−Δu = f in Ω

u = g on ∂Ω,

with Ω = (0, 1) × (0, 1) ⊂ R<sup>2</sup>. The Dirichlet data g and the right-hand side f are defined by

f(x, y) = 2π<sup>2</sup> cos(πx) cos(πy) and, g(x, y) = cos(πx) cos(πy).

**Functionality**

The main program is named as “***mgsolve.cpp***”. This file consists of 7 basic functionalities.

1) ***directsolve*** : This function calculates a single unknown **x** in a 3 by 3 matrix with the stencil and if matrix is not 3 by 3 , it performs Gauss-Seidel Solver and error calculation from using the functionality of the “***gssolve.hpp***” header file. 
1) ***Grid Prolongation***: This function extends the given grid **x** to a larger grid using a bi-linear interpolation prolongation technique, which is a fundamental operation in multigrid methods. Specifically, it takes a grid **x** and doubles its size along each dimension, effectively creating a finer grid. This finer grid is necessary for the multigrid V-Cycle method to iteratively refine the solution of PDEs.
1) ***Grid Restriction***: This function reduces the size of the given grid **x** to a coarser grid using a full weighting restriction technique. This operation is essential in multigrid methods like the V-Cycle method for iteratively solving PDEs.
1) ***Grid Residual***: This function calculates the residual of a given solution **x** with respect to the right-hand side **b** and the matrix **A** representing the coefficients of the PDE. This residual calculation is crucial for assessing the error in the current solution approximation. It iterates over each interior grid point in the provided solution **x**. Then it calculates the new approximation of the solution using neighboring points and stencil coefficients. Lastly, computation of the residual at each point by subtracting the right-hand side **b** from the newly calculated approximation is done.
1) ***Vcycle***: This function performs the steps for solving the PDE using restriction, prolongation, direct solve, and smoothing processes with Gauss-Seidel.
1) ***Errornorm***: This function performs discrete L2-norm calculation of the residual r = Au – b by following formula,

   ||r|| =  1n\*||r||<sub>2</sub>


1) ***main***: This function initializes the problem, sets the boundary conditions with parallization, calls the functions such as Vcycle, Residualnorm and saves the results. 

   In the first part, the argument-inserting process is done. It provides communication with the users in order to take the arguments. 

   In the second part, the initialization of the A, x, and b is done by the using Grid class and the definition of the boundary condition and inner points calculation is done.

   In the last part, Vcycle and Residual calculation functions are called to solve the problem. It should be noted that this part uses the *chrono* library for calculating the time that is spent by Vcycle and Residual functions. Also, this part saves the results in a solution file, additionally a report file. Lastly, it prints the runtime, convergence ratio, and L<sub>2</sub> norm to the console. 

**Dependencies**

- C++ compiler with OpenMP support (e.g. g++ compiler )

Example compilation command:

-O3 -fopenmp -Wall -Winline -Wshadow -std=c++17

- gssolve.hpp header file.

  This file includes functions such as constructing a grid class, point-wise summary function, Gauss-Seidel function, and residual norm calculation function.pp -fopenmp 


**Usage**

- **Command Line Arguments**: The program can be executed with command line arguments: number of levels and number of iterations.

Example: **./mgsolve 7 5** 

The first argument is the number of levels.

The second argument is the number of iterations.

If users do not provide a number of iteration, it is set automatically to 0. This means this V-Cycle program will iterate as much as possible until it converges a good result.

Example: **./mgsolve 7**

The first argument is the number of levels

The second argument is automatically set to 0 


- **Interactive Mode**: If no command line arguments are provided, the program prompts the user to input parameters interactively. This part is coded in the main function. In this mode, users are asked to insert a least 1 input, a number of levels. It is also asked to insert an iteration number optionally. 

Example: **./mgsolve** 

**>>This code needs at least 1 input argument to operate, the number of levels, would you like to insert it manually? Y/N**

**>> Y**

**>>Number of levels = *input(e.g. 4)***

**>> You also like to define the number of Vcycle iteration, put 0 to run until convergence = input (e.g. 2)**

The first argument is the number of levels

The second argument is automatically set to 0 if it is not inserted by the users.

**Outputs**

1) ***report.txt*** : Error, error norm, runtime values are logged in this file according to the number of level.
1) ***solution.txt*** : save the result of x, or unknown, according to the x and y location.

**Authors**

This program was authored by

Anish Adgaonkar @FAU

Arash Moradian @FAU

Mehmet Arif Bagci @ FAU

**Licence** 

This project is open-source. It is free to use, modify, and distribute this software.


