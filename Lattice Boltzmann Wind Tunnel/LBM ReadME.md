**Wind Tunnel Simulator using the Lattice Boltzmann Method**

**Overview**

This project implements a Lattice Boltzmann Method (LBM) for fluid dynamics simulations using a 2D grid of a wind tunnel. The code initializes the fluid domain, applies boundary conditions, computes collisions, and updates the lattice for each time step. The primary goal is to simulate fluid flow around obstacles and calculate drag and lift forces. The code has ability to read PGM file and compute the drag and lift of the obstacle in the wind tunnel. 

**Requirements**

- C++ compiler – g++ compiler
- MakeFile
- *parsing.hpp* header file
- OpenMP for parallel processing

**Functionality**

The main program *lbm.cpp* consists of 8  functions.

1) ***Domain*** : This function Filters the grid cells based on their flag (fluid, no-slip, inlet, outlet).
1) ***Pull*** : This function Updates fluid cells using neighboring cells' values.
1) **PullNoSlip**: This function updates boundary cells and computes drag and lift forces.
1) **VelocityCalculator**: This function computes velocities for fluid cells.
1) **InletInitializer**: This function sets the inlet boundary condition.
1) **OutletInitializer**: This function sets the outlet boundary condition.
1) **Collision**: This function computes the collision step of the LBM.
1) ***main*** : This function contains the logic for the Lattice Boltzmann simulation.


**Dependencies**

In order to use the program properly, the use *parsing.hpp* is mandatory.

- **Cell:** this structure represents an individual cell in the 2D grid with the following attributes:
1) flag: An integer flag indicating the cell's state (e.g., boundary conditions).
1) density: A double value representing the cell's density.
1) velocity: A vector of doubles representing the cell's velocity components.
1) boltz: A vector of doubles representing the Boltzmann distribution in various directions (C, N, E, S, W, NE, SE, SW, NW).
- **Grid**: this class manages a 2D array of Cell objects and provides various methods to manipulate and retrieve cell data.
- **write**: this function writes the grid data to a VTK file for visualization. This function outputs the grid's flags, density, and velocity for a specified timestep.
- **ObstacleFlag**: Flags cells within a radius r from the center (cx, cy) as obstacle cells (flag = 1)
- ` `**Grid PMGparsing**: Parses a PGM (Portable Gray Map) file to initialize the grid:

•  Reads the grid size from the PGM file.

• Flags cells based on the PGM values (cells with values less than 255 are flagged as obstacles).

- **Grid InputParsing**: Parses a parameter file to initialize the simulation
- Reads grid size, timesteps, velocity, Reynolds number, obstacle parameters, and VTK output settings from the file.
- Initializes the grid based on these parameters.
- Calculates the viscosity and relaxation time (tau) based on the input parameters.
- Flags obstacles if specified in the parameter file.

**Usage** 

The program accepts command-line arguments for setting parameters:

`	`**./lbm params.dat**

- params.dat= parameter file

If the user does not provide any file, the program shows message

`	`*Example*: **./lbm**

- Please provide the name of the parameter file for initialization

Accepted parameters in the parameter file are ;

- sizex: number of grid point in x direction in wind tunnel domain
- sizey: number of grid points in y direction in wind tunnel domain
- timesteps: the number of timesteps
- uin: the inlet velocity
- Re: Reynolds number used in calculating timesteps size
- spherex: x dimension of the center of the obstacle
- spherey: y dimension of the center of the obstacle
- vtk_file: name of the output file
- vtk_step: print interval of the output file 
- print_interval: print interval of the drag and lift (bonus task)
- geometry: input wind tunnel as a PGM file (bonus task)

**Output**

The output is vtk file based on the information that inserted in the parameter file.

Lift and Drag interval values will be printed in the screen with the interval inserted in the parameter file.  

**Licence**

This project is open-source. It is free to use, modify, and distribute this software.

**Author**

This program was authored by 

Anish Adgaonkar @ FAU

Arash Moradian @ FAU

Mehmet Arif Bagci @FAU

