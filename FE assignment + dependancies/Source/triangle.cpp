#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "Colsamm.h"
#include <map>
#include <math.h>

inline double k2(double x, double y){
return (100+0.01)*pow(M_E,(-50*(x*x+y*y))) -100; //sigma is set to 0.01
}

using namespace :: _COLSAMM_;

void Matrixcalculator(std::vector<std::map<int,double>> &Amatrix, std::vector<int> triangles, const std::vector<std::vector<double>> &vertices, bool mass)  //stiffness
 {
     std::vector<std::vector<double>> my_local_matrix;
std::vector<double> local_vetrices(6,0.0);
ELEMENTS::Triangle  my_element;
for (int i=0;i<vertices.size();i++){
    for (int j=0;j<vertices[i].size();j++){
        std::cout<<vertices[i][j]<< "  " ;
    }
    std::cout<<std::endl;
}

// array local_vertices contains the x and y coordinates of the
// triangle vertices in order, x0, y0, x1, y1, x2, y2
local_vetrices[0]=vertices[triangles[0]][0]; local_vetrices[1]=vertices[triangles[0]][1];
local_vetrices[2]=vertices[triangles[1]][0]; local_vetrices[3]=vertices[triangles[1]][1];
local_vetrices[4]=vertices[triangles[2]][0]; local_vetrices[5]=vertices[triangles[2]][1];

for (int i=0;i<local_vetrices.size();i++) {
    std::cout<<" "<<local_vetrices[i];
}
std::cout<<"\n local vertices ^" <<std::endl;
// pass the vertices to the finite element
my_element (local_vetrices);

// Select the equation to solve
if (mass=false) {
//calculates the local stiffness matrix
my_local_matrix=my_element.integrate(grad(v_())*grad(w_())-func<double>(k2)*v_()*w_());
}
else {
//calculates the local mass matrix
my_local_matrix=my_element.integrate(v_()*w_());
}

std::cout<<"local stiffness matrix" <<std::endl;
for( int i=0;i<my_local_matrix.size();i++){
    for (int j=0;j<my_local_matrix[i].size();j++){
        std::cout << i << " " << j << " " << my_local_matrix[i][j]<< " ";
    }
    std::cout<< std::endl;
}

//probably wrong
for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
        if (Amatrix[triangles[i]].count(triangles[j]) > 0) {
            Amatrix[triangles[i]][triangles[j]]+=my_local_matrix[i][j];
        }
        else {
            Amatrix[triangles[i]][triangles[j]]= my_local_matrix[i][j];        }
    }
}

std::cout<< "this is A matrix and this is the local matrix " << my_local_matrix[0][0]<<std::endl;


std::cout<< "I reached the end of function" <<std::endl;

}

int main(){
std::vector<std::vector<double>> vertices{{0,0},{0.25,0},{0,0.25},{0.25,0.25},{0,0.5},{0.25,0.5}};
std::vector<std::vector<int>> triangles={{0,1,2},{2,3,1},{2,3,5},{2,4,5}};
std::vector<std::map<int,double>> Amatrix(vertices.size()); //stiffness
std::vector<std::map<int,double>> Mmatrix; //mass
bool mass =false;

// from text file
for(unsigned long int i=0;i<Amatrix.size();i++){
    for (auto& j : Amatrix[i]){
        std::cout << i << " " << j.first << " " << j.second<< std::endl;
    }
}

std::cout<< "I reached the function call" <<std::endl;

for (int i=0;i<triangles.size();i++){
Matrixcalculator(Amatrix,triangles[i],vertices,false);
}
for(unsigned long int i=0;i<Amatrix.size();i++){
    for (auto& j : Amatrix[i]){
        std::cout << i << " " << j.first << " " << j.second;
    }
    std::cout<< std::endl;
}
std::cout<< "I reached the real end " <<std::endl;

return 0;
}
