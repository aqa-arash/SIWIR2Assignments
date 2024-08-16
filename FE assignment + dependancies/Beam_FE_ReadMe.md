**Beam Waveguide Propagation Solver**

**Overview**

This C++ program is designed to solve the propagation of information within a beam waveguide using finite element methods. It calculates the eigenmodes of the waveguide by solving a generalized eigenvalue problem. The code implements the inverse power iteration algorithm to find the eigenvalues and eigenvectors associated with the waveguide. This program reads the mesh from the "unit_circle.txt" file which should be in the same directory, to read other files you need to change the waveguide.cpp. Please ensure the formatting of the file is constant.

**Requirements**

- C++ compiler – optionally G++
- MakeFile
- Colsamm Library (to more details -> <https://www.cs10.tf.fau.de/research/software/expde/colsamm/> )

**Functionality** 

The main program *waveguide.cpp* consists of 5 basic functions. 

1) Double k2 : this function defines the gradient profile by the following equation. 

   k<sup>2</sup>(x, y) = (100 + δ)e−50(x<sup>2</sup>+y<sup>2</sup>) − 100.

1) ksq : this function prints and stores the k<sup>2</sup> values for each vertex in a file named *ksq.txt*.
1) mapfile : this function prints the sparse matrix stored as a map, the name of the file to store in is an input argument. Specifically this function is able to store the matrix for stiffness and mass.
1) Matrixcalculator: this function takes a sparse matrix stored as a map and populates it based on the given triangle, takes the vertices and a boolean value as input arguments. If the boolean value is true, it calculates the mass matrix. Otherwise, it calculates the stiffness matrix.
1) main : this function parses the input arguments to run the program. This also enables users to insert their specific values interactively. Firstly, it creates the vertices and triangles vectors and gets the *unit\_circle.txt* file to fetch the inputs.  Secondly, it creates the stiffness and mass matrix in a map format. It estimates the u<sub>h</sub> and normalizes these estimations with the Euclidean norm which is defined in *inversepower.hpp* file. It should also be noted that this function also enables to store of the u<sub>h</sub> values in the *eigenmode.txt* file. Thirdly, it populates the stiffness and mass matrix, for each triangle. Lastly, the lambda value is calculated by InversePower function which is defined in *inversepower.hpp* file. 

**Dependencies** 

In order to use the program properly, the use of *inversepower.hpp* and *io.hpp* files is mandatory.

1) *inversepower.hpp*: This code provides implementations of various numerical methods for solving systems of equations represented by sparse matrices, particularly focusing on matrices stored in map format. The provided functionalities include:

• Matrix-vector multiplication for sparse matrices. ( matrix\_vector\_multiplication function)  

• Calculation of the residual norm for a system of equations with a map coefficient matrix. (RNMap function) 

• Euclidean normalization of a given vector. (euclidean\_normalize function)

• Gauss-Seidel iterative method for solving systems of equations with a map coefficient matrix until a specified residual norm is achieved. (GSmap function)

• Dot product of two vectors. (vector\_vector\_mult function) 

• Inverse power algorithm for stiffness and mass matrices, utilizing Gauss-Seidel method as the solver. (Inversepower function)

1) io.hpp : this code reads the text files and stores the vertex and faces data in corresponding 2D vector. It only works with files that are formatted exactly as the example. In this case, *unit\_circle.txt*.

**Usage**

The program accepts command-line arguments for setting parameters:

- **sigma**: Fraction index related to the waveguide material.
- **epsilon**: Residual norm accepted for the Gauss-Seidel solver.
- **Writetofile**: Writes K<sup>2</sup>, Stiffness and Mass matrix into corresponding files if set to 1. By default it does not write those values to file.




**Command Line Arguments**: The program can be executed with command line arguments directly. 

Example: **./waveguide 0.01 0 1** 

- The first argument is the sigma, the second argument is the epsilon and the third argument is write to file status.

If the second argument is not inserted by the user, or it is equal to 0, it will automatically be set to 10<sup>-10</sup>. Please note if you want to write to file you need to provide a epsilon as well. 

Example: **./waveguide 0.01**

- The first argument is the the sigma second argument is set automatically to 10<sup>-10</sup>. By default it does not write K<sup>2</sup>, Stiffness and Mass values to file.

**Interactive Mode**: If no command line arguments are provided, the program prompts the user to input parameters interactively. This part is coded in the main function. In this mode, users are asked to insert a least 1 input, sigma. It is also asked to insert the convergence criteria, epsilon. Additionally, in this mode, it is possible to decide the third argument in the command line. In this part, users can decide to print ksq, stiffness, or mass matrix results to the file. 

Example: **./waveguide** 

**>>** **This code needs at least 1 input argument to operate, the value of sigma (related to the fraction index) = (0 to exit) = *input (e.g. 0.01)***

**>> You could also define the convergence criteria (residual norm limit) of GS (0 for default 10^-10) = *input (e.g. 0)***

**>> You could also write Stiffness and Mass matrix, as well as KSQ to file input 1 to write file = *input( e.g. A or M or KSQ)***

**Input Data**

The program runs with the input data in the form of a file named *unit\_circle.txt*, which includes information of vertices and triangles defining the waveguide geometry. In case of a change in the input, *waveguide.cpp* files should be modified. 

**Output**

- **eigenmode.txt**: Contains the eigenmodes of the waveguide.
- **ksq.txt**: (Optional) Contains the calculated k<sup>2</sup> values for each vertex.
- **A.txt** and **M.txt**: (Optional) Stiffness and mass matrices stored as sparse matrices.


**License**

This project is open-source. It is free to use, modify, and distribute this software.

**Authors**

This program was authored by 

Anish Adgaonkar @ FAU

Arash Moradian @ FAU

Mehmet Arif Bagci @FAU  

