#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <math.h>
#include <vector>
#include <string>
#include <sys/time.h>

#include "VTKparsing.hpp"

using namespace std;

//calculated the minimum and maximum x and y values for a set of vertices (bounding box)
std::vector<double> MinMax(std::vector<std::vector<double>> & vertices){
    double xmin=vertices[0][0],ymin=vertices[0][1],xmax=vertices[0][0],ymax=vertices[0][1];

    for (long unsigned int i=0;i<vertices.size();i++){
        if (xmin>vertices[i][0]) xmin=vertices[i][0];
        if (ymin>vertices[i][1]) ymin=vertices[i][1];
        if (xmax<vertices[i][0]) xmax=vertices[i][0];
        if (ymax<vertices[i][1]) ymax=vertices[i][1];
    }
    return {xmin,ymin,xmax,ymax};
}

//calculates barycentric coordinates of a node based on its global position within a triangle, given the index of vertices in the triangle and access to the vertices of the triangulation.
vector<double> BarycentricCoordCalc(const vector<int>& trianglesin,const vector<vector<double>>&vertices,const vector<double> &node) {
 // Extract vertex coordinates
    double x1 = vertices[trianglesin[0]][0], y1 = vertices[trianglesin[0]][1];
    double x2 = vertices[trianglesin[1]][0], y2 = vertices[trianglesin[1]][1];
    double x3 = vertices[trianglesin[2]][0], y3 = vertices[trianglesin[2]][1];
    double xp = node[0], yp = node[1];

    // Calculate the area of the entire triangle using determinant method
    double area = 0.5 * fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

    // Calculate the area of sub-triangles formed with the node
    double area1 = 0.5 * fabs((x2 - xp) * (y3 - yp) - (x3 - xp) * (y2 - yp));
    double area2 = 0.5 * fabs((x1 - xp) * (y3 - yp) - (x3 - xp) * (y1 - yp));
    double area3 = 0.5 * fabs((x1 - xp) * (y2 - yp) - (x2 - xp) * (y1 - yp));

    // Calculate barycentric coordinates as ratios of sub-triangle areas to the total area
    double u = area1 / area;
    double v = area2 / area;
    double w = area3 / area;

    return {u, v, w};
}

// given a triangulation and its underlying vertices, as well as the chosen size of the auxiliary grid and its bounding box,
// populates the auxiliary grid with corresponding triangles
vector<vector<int>> PopulateStructuredGrid(const vector<vector<int>>&trianglesin,
                                           const vector<vector<double>> &verticesin,
                                           double h,const vector<double> &minmaxcoords){
//calculate number of nodes in x and y direction in the auxiliary grid
    int nx = 1+(minmaxcoords[2]-minmaxcoords[0])/h;
    int ny = 1+(minmaxcoords[3]-minmaxcoords[1])/h;

//declare a 1D vector corresponding to each cell in the auxiliary grid
    vector<vector<int>> structuredcells(nx*ny);

//iterate through triangles in the source triangulation
    for (long unsigned int k=0;k<trianglesin.size();k++){
        vector<double> localminmax;
        vector<vector<double>> trianglenodes(3);


//create the set of vertices of the triangle to pass to the MinMax Function
        trianglenodes[0].push_back(verticesin[trianglesin[k][0]][0]);
        trianglenodes[0].push_back(verticesin[trianglesin[k][0]][1]);
        trianglenodes[1].push_back(verticesin[trianglesin[k][1]][0]);
        trianglenodes[1].push_back(verticesin[trianglesin[k][1]][1]);
        trianglenodes[2].push_back(verticesin[trianglesin[k][2]][0]);
        trianglenodes[2].push_back(verticesin[trianglesin[k][2]][1]);

//calculate the bounding box of the triangle
        localminmax=MinMax(trianglenodes);

//calculate the span of the bounding box
        int jfirst=(localminmax[1]-minmaxcoords[1])/h;
        int jmax=(localminmax[3]-minmaxcoords[1])/h;
        int ifirst=(localminmax[0]-minmaxcoords[0])/h;
        int imax=(localminmax[2]-minmaxcoords[0])/h;

//iterate through the span of the bounding box, populating the corresponding cells with the index of the triangle
        for (int j=jfirst;j<jmax+1;j++){
            for (int i=ifirst;i<imax+1;i++){
                structuredcells[(j*nx)+i].push_back(k);
                }
            }
        }
    return structuredcells;
}

//program is designed to take in the name of the source and target unstructured grids and interpolate between them using an auxiliary structured grid
 int main (int argc,char** argv){

    std::string targetfile = "./sample_grids/mesh_shifted_out.vtk";
    std::string infile = "./sample_grids/mesh_shifted_in.vtk";
    std::string outfile = "./sample_grids/mesh_shifted_out_interpolated.vtk";
    double h=0;
    bool print=false;
    if (argc==4){
        infile=argv[1];
        targetfile=argv[2];
        outfile = argv[2];
        outfile.insert(outfile.length()-4,"_interpolated");
        //cout<<"Result of the interpolation will be saved in " <<outfile<<endl;
        int power=stoi(argv[3]);
        h=1/pow(2,power);
        print=true;
        }
    else if (argc==3){
        infile=argv[1];
        targetfile=argv[2];
        outfile = argv[2];
        outfile.insert(outfile.length()-4,"_interpolated");
        //cout<<"Result of the interpolation will be saved in " <<outfile<<endl;
        }
    else {
        cout<<"The program can run with 'source file name' 'target file name' and optional h (auxiliary grid size) commandline arguments" <<endl;
        cout<< "Running interpolation between " <<infile<< " and "<< targetfile<<endl;
        cout<<"Result of the interpolation will be saved in " <<outfile<<endl;
        }


//containers
    std::vector<double> MinMaxCoords(4),temp(4);
    std::vector<std::vector<double>> verticesin,verticesout;
    std::vector<std::vector<double>> valuesin,valuesout;
    std::vector<std::vector<int>> trianglesin,trianglesout;

//read input and target files
    read(verticesin,trianglesin,valuesin,infile);
    read(verticesout,trianglesout,valuesout,targetfile);

cout<< "Input read #nodes " << verticesin.size()<<" #values "<< valuesin.size()<<" #cells "<< trianglesin.size()<< endl;
cout<< "Output read #nodes " << verticesout.size()<<" #values "<< valuesout.size()<<" #cells "<< trianglesout.size()<< endl;

//default value of h, if not declared by the user
    if (h==0) h=1.0/(sqrt(verticesin.size()));

//calculating the bounding box of the in and target grids
    MinMaxCoords=MinMax(verticesin);
    temp=MinMax(verticesout);

//choosing the bounding box of the auxiliary grid so that it encompasses both in and target grid
    if (MinMaxCoords[0]>temp[0]) MinMaxCoords[0]=temp[0];
    if (MinMaxCoords[1]>temp[1]) MinMaxCoords[1]=temp[1];
    if (MinMaxCoords[2]<temp[2]) MinMaxCoords[2]=temp[2];
    if (MinMaxCoords[3]<temp[3]) MinMaxCoords[3]=temp[3];

//number of nodes in the auxiliary grid in x and y direction
    int structuredgridpointsx = 1+(MinMaxCoords[2]-MinMaxCoords[0])/h;
    //int structuredgridpointsy = 1+(MinMaxCoords[3]-MinMaxCoords[1])/h;


//cout<< "Building structured grid with size "<< structuredgridpointsx*structuredgridpointsy<<endl;
//cout << "xmin " << MinMaxCoords[0]<<" xmax "<< MinMaxCoords[2] <<" ymin " << MinMaxCoords[1] <<" ymax "<< MinMaxCoords[3]<<endl;
    struct timeval t0, t;
    gettimeofday(&t0, NULL);
    vector<vector<int>> structuredcells=PopulateStructuredGrid(trianglesin,verticesin,h,MinMaxCoords);
//cout<< "Structured grid populated "<< endl;

//cout<< "start interpolation "<<endl;
//interpolation, going through each vertex in target grid
    for (long unsigned int i=0;i<verticesout.size();i++){
        vector<double> node= {verticesout[i][0],verticesout[i][1]};
    //calculate the location of the node in the auxiliary grid
        int xcoord= (node[0]-MinMaxCoords[0])/h;
        int ycoord= (node[1]-MinMaxCoords[1])/h;

        long unsigned int k=0;
        vector<double> BarycentricCoords(3);
        //if (i%10==0) cout<<'*';
        //iterate through triangles in auxiliary grid
            for (;k<structuredcells[xcoord+ycoord*structuredgridpointsx].size();k++){
            BarycentricCoords= BarycentricCoordCalc(trianglesin[structuredcells[xcoord+ycoord*structuredgridpointsx][k]],verticesin,node);
            //check if the point is inside the triangle
            if (((BarycentricCoords[0]+BarycentricCoords[1]+BarycentricCoords[2]-1)<0.000000001) &&
                BarycentricCoords[0] >0 && BarycentricCoords[0]<1 &&
                BarycentricCoords[1] >0 && BarycentricCoords[1]<1 &&
                BarycentricCoords[2] >0 && BarycentricCoords[2]<1 ) break;
            else if (k==structuredcells[xcoord+ycoord*structuredgridpointsx].size()-1) cout<<"point "<< i<<" was not in any source triangles!"<<endl;
            }
    //interpolation over the triangle using the barycentric coordinates
        valuesout[i][0]=   valuesin[trianglesin[structuredcells[xcoord+ycoord*structuredgridpointsx][k]][0]][0] * BarycentricCoords[0] +
                        valuesin[trianglesin[structuredcells[xcoord+ycoord*structuredgridpointsx][k]][1]][0] * BarycentricCoords[1] +
                        valuesin[trianglesin[structuredcells[xcoord+ycoord*structuredgridpointsx][k]][2]][0] * BarycentricCoords[2] ;
}

    gettimeofday(&t, NULL);
    std::cout << "\nWall clock time of interpolation: " <<
    ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t)1000000 +
    (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3
    << " ms" << std::endl;

 if (print==true){
        // To facilitate testing, outputs test results to file.
    ofstream report;
    report.open ("error.txt",std::ios_base::app);
        //report << setw(20)<< "run time" << setw(20)<< "Error norm"<< setw(20)<<"gridsize"<<setw(20)<<"number of cores"<<endl;
    report <<setw(20)<<setprecision(15)<<h<<setw(20)<<setprecision(15)<< ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t)1000000 +
                                                                                (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3
                                                                                << " ms"  <<endl;
    report.close();
    }
//write out file
    write(verticesout,trianglesout,valuesout,outfile);
    return 0;
 }



