**VTK Grid Interpolation**

**Overview**

This C++ program is designed to perform interpolation between two unstructured grids, stored in VTK file format, using an auxiliary structured grid. The goal is to map values from a source grid to a target grid. The program reads the grids, constructs an auxiliary grid, computes barycentric coordinates, and interpolates values for the target grid nodes based on the source grid triangles. Thanks to this program, it is possible to carries out an optimized interpolation of scalar values from one unstructured grid to another.

**Requirements**

- C++ compiler – optionally G++
- MakeFile
- _VTKparsing.hpp_ header file for reading from VTK file and writing a VTK file.

**Functionality**

The main program _intug.cpp_ consists of 4 basic functions.

1. **_MinMax_** : This function calculates the minimum and maximum x and y values for a given set of vertices, effectively determining the bounding box for the vertices. This function gets a vector of vertex coordinates as a parameter, where each vertex is represented as a vector of doubles \[x,y\].
2. **_BarycentricCoordCalc_** : This function calculates the barycentric coordinates for a given node within a triangle, using the triangle’s vertex indices and coordinates.
3. **_PopulateStructuredGrid_** : This function creates and populates an auxiliary structured grid with triangles from the source grid, aiding in efficient interpolation.
4. **_main_** : The main function controls the flow of the program, managing input/output operations, grid construction, and interpolation.

First, The program reads the source and target files. It also computes the auxiliary grid size h if provided or uses a default value. then, it reads vertices, triangles, and values from the source and target VTK files using the read function provided in VTKparsing.hpp.
afterwards, it computes the bounding boxes for the source grids and determines a combined bounding box for the auxiliary grid.  The auxiliary grid is, then, constructed using PopulateStructuredGrid, which associates source Ingrid triangles with the auxiliary grid cells.
Finally, the program performs interpolation for each vertex in the target grid, computing the barycentric coordinates and interpolating values based on the triangles in the corresponding auxiliary grid cell and interpolated Grid is written to the output file using the write function from VTKparsing.hpp.

*please note that the OUTGRID should be completely within the INGRID, as the program is designed to interpolate and not to extrapolate*
**_For the bonus task_**, if the optional parameter h is provided, the runtime is logged to an error.txt file for performance analysis.

**Dependencies**

In order to use the program properly, the use _VTKparsing.hpp_ is mandatory.

1. **_VTKparsing.hpp_** : This header file contains two C++ functions to handle the reading and writing of VTK file formats. These functions are;;

• _read_: this function reads a VTK file and extracts the vertices, faces, and point data (scalar values). It interprets the VTK file by checking keywords such as POINTS, CELLS, and POINT_DATA to determine which type of data it is reading.

The function reads the file line by line:

For POINTS, it reads the vertex coordinates.

For CELLS, it reads the vertex indices that form each triangle.

For POINT_DATA, it reads the scalar values associated with each vertex.

• _write_: this function creates a VTK file from the given vertices, faces, and point data. It writes the data in a format that can be used for visualization in VTK-supported tools.

**Usage**

The program accepts command-line arguments for setting parameters:

- INGRID.vtk: the grid containing scalar values on each grid node.
- OUTGRID.vtk: the grid to which the programmers should interpolate these scalar values.
- h= mesh size (optional)

**Command Line Arguments**: The program can be executed with command line arguments directly.

Example: **./intug mesh_fine_in.vtk mesh_fine_out.vtk**

- The first argument is the ingrid and the second argument is the outgrid. If h size is not inserted, it is automatically set to 1/sqrt(N) where N is the number of nodes in the INGRID.

**Interactive Mode**: If no command line arguments are provided, the program runs with the hardcoded INGRID and OUTGRID files. The output of this mode as follows.  

Example: **./intug**

**\>>** **This program can run with ‘source file name’ ‘target file name’ and optional h (auxiliary grid size) command line arguments.**

**Running interpolation between ./sample_grids/mesh_shifted_in.vtk and ./sample_grids/mesh_shifted_out.vtk**

**Result of the interpolation will be saved in ./sample_grids/mesh_shifted_out_interpolated.vtk**

**Input Data**

The program runs with the inputs data from the directory. In case of failure of the program, please be sure all the files are in the same directory with the program and check the path from the code and adjust it according to your file order.

**Output**

The program gives outputs in the following format: **_OUTGRID_interpolated.vtk_**

**License**

This project is open-source. It is free to use, modify, and distribute this software.

**Authors**

This program was authored by

Anish Adgaonkar @ FAU

Arash Moradian @ FAU

Mehmet Arif Bagci @FAU