#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iomanip>
using namespace std;
// Grid class to control 2D arrays
class Grid
{
public:
  // allocate memory in constructor
  Grid(int __nx__, int __ny__)
      : nx(__nx__)
      , ny(__ny__)
  {
      data = std::vector<double>(nx * ny); // Initialize vector with size
  }

  // free memory in destructor (not required for std::vector)
  ~Grid() {} // Empty destructor

  // Quick assign method
    void quickset(std::vector<vector<double>> input)
    {   int k=0;
        for (int j=0;j<ny;j++) {
            for (int i=0;i<nx;i++) {
                data[k]=input[i][j];
                k++;
            }
        }
    }

    //even quicker - unity
     void evenquickerset(double input)
    {
        #pragma omp parallel for
        for (int k=0;k<ny*nx;k++) {
           data[k]=input;
        }
    }

    // to quickly fill the boundary with a value
    void boundaryset (double x){
    for (int i=0;i<nx;i++){
        data[i + ny * nx]=x;
        data[i]=x;}
    for (int j=0;j<ny;j++){
        data[j*nx]=x;
        data[(j+1)*nx]=x;
    }
    }


  // getter
  double operator()(int i, int j) const { return data[i + j * nx]; }

  // setter
  double& operator()(int i, int j) { return data[i + j * nx]; }

  // access dimensions- Domain size is considered 1x1
  int nyp() const { return ny; }//number of points in y direction
  int nxp() const { return nx; }// number of points in x direction
  double hx() const {return 1.0/(nx-1);}
  double hy() const {return 1.0/(ny-1);}


private:
  // internal data
  int nx, ny;
  std::vector<double> data;
};

//adds the elemnts of B, to corresponding elements in x 
void Pointwisesum (Grid&A,const Grid& B){
    for (int j=0;j<A.nyp();j++){
        for (int i=0;i<A.nxp();i++){
            A(i,j)=A(i,j)+B(i,j);
        }
    }
}

//Gauss seidel function, changes the value of initial estimate x, based on b,A and the number of iterations. where A is a stencil
void GSS(Grid& x,const Grid& b,int iter,const Grid& A){
//solve Ax=b with Gauss Seidel method

double C,W,E,N,S,NW,NE,SW,SE;
// calculating the stencil
SW=A(0,0)*(1.0/(x.hx()*x.hy()));
S= A(1,0)*(1.0/pow(x.hy(),2));
SE=A(2,0)*(1.0/(x.hx()*x.hy()));
W=A(0,1)*(1.0/pow(x.hx(),2));
E=A(2,1)*(1.0/pow(x.hx(),2));
NW=A(0,2)*(1.0/(x.hx()*x.hy()));
N=A(1,2)*(1.0/pow(x.hy(),2));
NE=A(2,2)*(1.0/(x.hx()*x.hy()));
C=(A(1,1)*((1.0/pow(x.hx(),2))+(1.0/pow(x.hy(),2))));


// pow for square power

for (int k=0;k<iter;k++){ //main iteration
    //------red
    int i,j;
    #pragma omp parallel for private(i)
    for ( j=1; j<x.nyp()-1;j++)// y direction iteration
        {
            if(j%2==0){
                   i=1; //jumping over
            }
            else {
                i=2; // switching the jump
            }

			for (;i<x.nxp()-1;i+=2){//iterate x direction
				double newparts=0.0; // repository for the sum values
				newparts+=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1); //calculate based on old estimates
				x(i,j)=((b(i,j)- newparts) /C); //calculate new estimate of xi
			}
        }

    //-----black
    #pragma omp parallel for private(i)
    for (j=1; j<x.nyp()-1;j++)// y direction iteration
        {
            if(j%2!=0){
                   i=1; //jumping 1 point
            }
            else {
                i=2; // switching on the next line
            }

			for (;i<x.nxp()-1;i+=2){//iterate x direction
				double newparts=0.0; // repository for sum
				newparts+=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1); //calculate based on old estimates
				x(i,j)=((b(i,j)- newparts) /C); //calculate new estimate of xi
			}
        }
}
}

// calculates the norm of the residual error, based on x,b and A
double Residualnorm (const Grid& x,const Grid& b,const Grid& A){
double C,W,E,N,S,NW,NE,SW,SE;
double er=0.0;
double n=0;//number of grid nodes
// calculating the stencil
SW=A(0,0)*(1.0/(x.hx()*x.hy()));
S= A(1,0)*(1.0/pow(x.hy(),2));
SE=A(2,0)*(1.0/(x.hx()*x.hy()));
W=A(0,1)*(1.0/pow(x.hx(),2));
E=A(2,1)*(1.0/pow(x.hx(),2));
NW=A(0,2)*(1.0/(x.hx()*x.hy()));
N=A(1,2)*(1.0/pow(x.hy(),2));
NE=A(2,2)*(1.0/(x.hx()*x.hy()));
C=(A(1,1)*((1.0/pow(x.hx(),2))+(1.0/pow(x.hy(),2))));

#pragma omp parallel for reduction(+:er)
for (int j=1; j<x.nyp()-1;j++)// y direction iteration
    {
        for (int i =1 ;i<x.nxp()-1;i++){//iterate x direction
            double newparts=0.0;// repository for the Ax
            newparts=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1)+C*x(i,j); //calculate based on old estimates
            newparts=(b(i,j)-newparts); //adds the second degree of error term for each point
			n++;
            er+=newparts*newparts;//r square sum
            }
    }
//returns square root of the sum above, divided by the number of grid points.
return sqrt(er/n);
}

//writes the solution to the text file
void Writesolution(const Grid &x){
    ofstream solutionfile;
    solutionfile.open ("solution.txt");
    for (int j=x.nyp()-1; j>=0 ;j--) { //printing results of x
        for (int i=0;i<x.nxp();i++){
               solutionfile << i*x.hx() <<" " <<j*x.hy() << " "<<x(i,j)<<"\n";
        }
    }
    solutionfile.close();
cout<<"The result is available in the file solution.txt \n";	
}