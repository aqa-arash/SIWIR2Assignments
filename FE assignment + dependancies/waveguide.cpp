#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <vector>
#include <map>
#include <math.h>
#include "./Source/Colsamm.h"
#include "io.hpp"
#include "inversepower.hpp"

using namespace :: _COLSAMM_;

//global variable used in k^2 calculation
double sigma=0.01;

// k^2 function
double k2(double x, double y){
return (100+sigma)*pow(M_E,(-50*(x*x+y*y))) -100; //sigma is set to 0.01
}

//prints the K^2 value for each vertex to a file named ksq.txt
void ksq(std::vector<std::vector<double>> &vertices){
    ofstream textfile;
    textfile.open ("ksq.txt");
    for (long unsigned int i=0;i<vertices.size();i++){
        textfile <<vertices[i][0]<<" "<<vertices[i][1]<<" "<< k2(vertices[i][0],vertices[i][1]) <<endl;
    }
    textfile.close();
}

//prints the sparse matrix stored as map, name of the file to store in is an input argument
void mapfile(std::vector<std::map<int,double>> &matrix,std::string a){
    ofstream textfile;
    std::string b= a+".txt";
    textfile.open (b);
    for (long unsigned int i=0;i<matrix.size();i++){
        for (auto& j : matrix[i]){
            textfile <<i<<" "<<j.first<<" "<< j.second <<endl;
            }
        }
    textfile.close();
}

//takes a sparse matrix stored as map and populates it based on the given triangle, takes the vertices and a bool value as input arguments
//if bool is true, calculates mass matrix, else calculates stiffness
void Matrixcalculator(std::vector<std::map<int,double>> &matrix, std::vector<int> triangles, const std::vector<std::vector<double>> &vertices, bool mass){
    std::vector<std::vector<double>> my_local_matrix;
    std::vector<double> local_vetrices(6,0.0);
    ELEMENTS::Triangle  my_element;


    // array local_vertices contains the x and y coordinates of the
    // triangle vertices in order, x0, y0, x1, y1, x2, y2
    local_vetrices[0]=vertices[triangles[0]][0]; local_vetrices[1]=vertices[triangles[0]][1];
    local_vetrices[2]=vertices[triangles[1]][0]; local_vetrices[3]=vertices[triangles[1]][1];
    local_vetrices[4]=vertices[triangles[2]][0]; local_vetrices[5]=vertices[triangles[2]][1];

    // pass the vertices to the finite element
    my_element (local_vetrices);

    // Select the equation to solve
    if (mass==false) {
        //calculates the local stiffness matrix
        my_local_matrix=my_element.integrate(grad(v_())*grad(w_())-func<double>(k2)*v_()*w_());
        }
    else {
        //calculates the local mass matrix
        my_local_matrix=my_element.integrate(v_()*w_());
        }

    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            if (matrix[triangles[i]].count(triangles[j]) > 0) {
                matrix[triangles[i]][triangles[j]]+=my_local_matrix[i][j];
                }
            else {
                matrix[triangles[i]][triangles[j]]= my_local_matrix[i][j];
                }
            }
        }
}

//solve propagation of information withing a beam waveguide
int main(int argc,char**argv){
double epsilon;
int writetofile=0;
// sigma is the fraction index and epsilon is residual norm accepted for GS
//if epsilon = 0, runs until residual norm is 10^-10
//if write to file ==1 => writes KSQ, A and M matrix
//input argument parsing
if (argc==4){
 sigma=stod(argv[1]);
 epsilon=stod(argv[2]);
  if (epsilon==0){
    epsilon=0.0000000001;
 }
 writetofile=stoi(argv[3]);
}
else if (argc==3){
 sigma=stod(argv[1]);
 epsilon=stod(argv[2]);
 if (epsilon==0){
    epsilon=0.0000000001;
 }
}
else if (argc==2){
     sigma = stod(argv[1]);
     std::cout<<"No convergence criteria was defined for the GS solver, 10^-10 is selected by default \n";
     epsilon =0.0000000001;
}
else {
    cout<< "This code needs at least 1 input argument to operate, the value of sigma (related to the fraction index) = (0 to exit) \n ";
    cin>>sigma;
    if (cin.fail()) {
          cout << "That was not a double." << endl;
          return -1;
                }
    if (sigma==0){
        std::cout<< "See you next time! \n";
        return -1;
    }
    cout<< "You could also define the convergence criteria (residual norm limit) of GS (0 for default 10^-10) = ";
    cin>>epsilon;
    if (cin.fail()) {
        cout << "That was not a double." << endl;
        return -1;
        }
    if (epsilon==0){
        epsilon =0.0000000001;
        }
    cout<< "You could also write Stiffness and Mass matrix, as well as KSQ to file input 1 to write file = ";
    cin>>writetofile;
      if (writetofile!=1){
        writetofile=0;
        }
    }
// input argument parsing done


std::vector<std::vector<double>> vertices;
std::vector<std::vector<int>> triangles;
input(vertices,triangles,"unit_circle.txt");//reads the information from the file "unit_circle.txt" in the same directory and stores it in the vectors above



std::vector<std::map<int,double>> Amatrix(vertices.size()); //stiffness
std::vector<std::map<int,double>> Mmatrix(vertices.size()); //mass
std::vector<double> uh(vertices.size()); //resulting vector of beam propagation

//initial estimate of uh all 1
#pragma omp parallel for
for (long unsigned int i =0; i<uh.size();i++){
    uh[i]= 1;
}
euclidean_normalize(uh); //normalize the initial estimate

// populates the stiffness and mass matrix, for each triangle
for (long unsigned int i=0;i<triangles.size();i++){
    Matrixcalculator(Amatrix,triangles[i],vertices,false);
    Matrixcalculator(Mmatrix,triangles[i],vertices,true);
}

// write the stiffness and mass matrix to files A.txt and M.txt and K^2 to ksq.txt,
if (writetofile==1){
ksq(vertices);
mapfile(Amatrix,"A");
mapfile(Mmatrix,"M");
}

//running the inverse power algorithm from hpp file
double lamda=Inversepower(Amatrix,Mmatrix,uh,epsilon);

//printing uh values in the eigenmode.txt file
ofstream textfile;
textfile.open ("eigenmode.txt");
for (long unsigned int i=0;i<vertices.size();i++){
     textfile <<vertices[i][0]<<" "<<vertices[i][1]<<" "<< uh[i] <<endl;
    }
textfile.close();


std::cout<< "The final lamda value is "<< lamda <<std::endl;

return 0;
}
