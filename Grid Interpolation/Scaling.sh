#!/bin/bash
echo 
echo "testing grid interpolation with various auxiliary grid size"
echo "Auxiliary grid size is 1/(2^h)"
echo "results will be written to the report.txt file"
./intug mesh_finest_in.vtk mesh_finest_out.vtk 0
echo h=0 done
./intug mesh_finest_in.vtk mesh_finest_out.vtk 1
echo h=1 done
./intug mesh_finest_in.vtk mesh_finest_out.vtk 2
echo h=2 done
./intug mesh_finest_in.vtk mesh_finest_out.vtk 3
echo h=3 done
./intug mesh_finest_in.vtk mesh_finest_out.vtk 4
echo h=4 done
./intug mesh_finest_in.vtk mesh_finest_out.vtk 5
echo h=5 done
./intug mesh_finest_in.vtk mesh_finest_out.vtk 6
echo h=6 done
./intug mesh_finest_in.vtk mesh_finest_out.vtk 7
echo h=7 done
./intug mesh_finest_in.vtk mesh_finest_out.vtk 8
echo h=8 done
./intug mesh_finest_in.vtk mesh_finest_out.vtk 9
echo h=9 done
echo bash file done