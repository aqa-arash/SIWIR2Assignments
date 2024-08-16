#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstring>
#include <iomanip>
#include "gssolve.hpp"
using namespace std;

// tries to solve the 3x3 system of eqations directly (if handed in a non 3x3 system, runs GS until convergence
void directsolve (Grid& x,const Grid& b,const Grid& A){
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
//calculate the single unknown (middle node)
if (x.nyp()==x.nxp() && x.nyp()==3) {
    double newparts=0.0;
    newparts+=SW*x(0,0)+S*x(1,0)+SE*x(2,0)+W*x(0,1)+E*x(2,1)+NW*x(0,2)+N*x(1,2)+NE*x(2,2); //calculate based on old estimates
    x(1,1)=((b(1,1)- newparts) /C);
}
else { //if x is not 3x3, runs gs until convergence
    double er=1.0;
	while (er>pow(10,-9)){
    GSS(x,b,5,A);
    er=Residualnorm(x,b,A);
}
}}


//prolongation to a finer grid using bilinear interpolation 
Grid Prolongation(const Grid& x){

int nx=2*x.nxp()-1;
int ny=2*x.nyp()-1;
Grid finer(nx,ny);
finer.evenquickerset(0);
#pragma omp parallel for
for (int j=2;j<ny;j+=2){
    for (int i=2;i<nx;i+=2){
        finer(i,j)=0.5*x(i/2,j/2);
    }
    for (int i=1;i<nx;i+=2){
        finer(i,j)=0.25*(x(i/2,j/2)+x((i/2)+1,j/2));
    }
}
#pragma omp parallel for
for (int j=1;j<ny;j+=2){
    for (int i=2;i<nx;i+=2){
        finer(i,j)=0.25*(x(i/2,j/2)+x((i/2),(j/2)+1));
    }
    for (int i=1;i<nx;i+=2){
        finer(i,j)= 0.125*(x(i/2,j/2)+x(i/2,(j/2)+1)+x((i/2)+1,(j/2)+1)+x((i/2)+1,(j/2)));
    }
} 
return finer;
}

//restriction to a coarser grid using full weighting
Grid  Restriction(const Grid& x){

int nx=x.nxp()/2;
int ny=x.nyp()/2;
Grid coarser(nx+1,ny+1);
coarser.evenquickerset(0);
#pragma omp parallel for
for (int j=1;j<ny;j++){
    for (int i=1;i<nx;i++){
        double sumcardinal,sumordinal;
        sumcardinal= 0.25*(x(2*i-1,2*j)+x(2*i+1,2*j)+x(2*i,2*j+1)+x(2*i,2*j-1));
        sumordinal=0.125*(x(2*i+1,2*j+1)+x(2*i+1,2*j-1)+x(2*i-1,2*j+1)+x(2*i-1,2*j-1));
        coarser(i,j)=0.5*x(2*i,2*j)+sumcardinal+sumordinal;
    }
}
return coarser;
}


// calculates the residual error, based on x,b and A
Grid Residual (const Grid& x,const Grid& b,const Grid& A){
double C,W,E,N,S,NW,NE,SW,SE;
Grid er(x.nxp(),x.nyp());
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

er.evenquickerset(0); //sets the er for all nodes to zero
#pragma omp parallel for
        for (int j=1; j<x.nyp()-1;j++)// y direction iteration
        {
			for (int i =1 ;i<x.nxp()-1;i++){//iterate x direction
				double newparts=0.0;// repository for the Ax
				newparts=SW*x(i-1,j-1)+S*x(i,j-1)+SE*x(i+1,j-1)+W*x(i-1,j)+E*x(i+1,j)+NW*x(i-1,j+1)+N*x(i,j+1)+NE*x(i+1,j+1)+C*x(i,j); //calculate based on old estimates
				er(i,j)=(b(i,j)-(newparts)); //adds the second degree of error term for each point
            }
        }
return er;}

// V cycle method, using the GS(2,1) as smoother, lvl(first input) is the number of levels and needs to correspond to the size of x and b
void Vcycle (int lvl,Grid& x,const Grid& b, const Grid& A){
	
std::vector<Grid> X,B; //creates repository for x and b on each level
for (int i=0;i<lvl+1;i++){
    if (i==0) { //first step, the actual Ax=b
        X.push_back(x);
        B.push_back(b);
    }
    else {
        X.push_back(Grid(pow(2,lvl-i+1)+1,pow(2,lvl-i+1)+1)); //define the coarser grid
        X[i].evenquickerset(0);//new left hand side estimate
        B.push_back(Restriction(Residual(X[i-1],B[i-1],A))); //populate the new right hand side with residuals of the former step
    }

    if (i==lvl){
		directsolve(X[i],B[i],A); //solve the 3x3 directly
    }
    else {
        GSS(X[i],B[i],2,A);// smoothing
    }
}
for (int i=lvl;i>0;i--){ //going up the Vcycle
    Pointwisesum(X[i-1],Prolongation(X[i]));//prolongating the results of the coarser grid and adding it to the finer
    GSS(X[i-1],B[i-1],1,A);//smoothing
}
x=X[0];//assign the final estimation to x
}

// calculates error norm for this particular problem, when given x as estimate
double Errornorm (const Grid& x){
double errornorm=0,pi=3.141592653589793238462643383279502884;

#pragma omp parallel for reduction (+:errornorm)
for (int j=0;j<x.nyp();j++){
    for (int i=0;i<x.nxp();i++){
        double error=0.0;
        error = cos(pi*i*x.hx())*cos(pi*j*x.hy()) - x(i,j);
        errornorm+=error*error;
    }
}
return sqrt(errornorm/(x.nyp()*x.nxp()));
}


int main(int argc,char** argv){

int lvl,iter;
// lvl is the number of levels and iter is iterations of Vcycle(2,1)  lvl=0-> 3x3 initial matrix
//if iter = 0, runs until residual norm is 10^-11
if (argc==3){
 lvl=stoi(argv[1]);
 iter=stoi(argv[2]);
}
else if (argc==2){
     lvl = stoi(argv[1]);
     iter =0;
}
else {
    cout<< "This code needs at least 1 input argument to operate, the number of levels, would you like to insert it manually? Y/N \n ";
    char A;
    cin >> A;
    if (A=='Y' || A=='y'){
          cout << "\n Number of levels = ";
          cin>>lvl;
          if (cin.fail()) {
                    cout << "That was not an integer." << endl;
                    return -1;
                }
        cout<< "You also like to define the number of Vcycle iteration, put 0 to run until convergence = ";
        cin>>iter;
          if (cin.fail()) {
                    cout << "That was not an integer." << endl;
                    return -1;
                }}

}

//variable definition
Grid stencil(3,3); //define variables as Grid A for the stencil
Grid x(pow(2,lvl+1)+1,pow(2,lvl+1)+1);
Grid b(pow(2,lvl+1)+1,pow(2,lvl+1)+1);
double er; // for error
stencil.quickset({{0,-1,0},{-1,2,-1},{0,-1,0}});// Coefficient matrix
x.evenquickerset(0);//x vector (initial estimate)
b.evenquickerset(0);//b vector (right hand side)
double pi=3.141592653589793238462643383279502884;

//define the boundry condition - basically g function
//in x direction
#pragma omp parallel for
for(int i=0;i<x.nxp();i++ ){
    x(i,x.nyp()-1)=-cos(pi*i*x.hx());
    x(i,0)=cos(pi*i*x.hx());
}
#pragma omp parallel for
// in y direction
for(int i=0;i<x.nyp();i++){
    x(0,i)=cos(pi*i*x.hy());
    x(x.nxp()-1,i)=-cos(pi*i*x.hy());
}

//define the b for the internal point -- basically f function
#pragma omp parallel for
for(int i=0;i<x.nxp();i++ ){
        for(int j=0; j<x.nyp();j++){
            b(i,j)= 2*(pi*pi)*cos(pi*i*x.hx())*cos(pi*j*x.hy());
        }
}

//calculate initial residual norm
double res1=Residualnorm(x,b,stencil);
double res2=0.0;

 //start of the V cycle process
if (iter!=0){
	auto start2= std::chrono::system_clock::now(); //wall clock
	for (int k=0;k<iter;k++){
		Vcycle(lvl,x,b,stencil);
		res2=Residualnorm(x,b,stencil);
		cout<<"New L2 Norm of residual = " <<res2<<" The Convergance ratio for this iteration  "<< res2/res1 <<endl;
		res1=res2;
	}
	auto fin2= std::chrono::system_clock::now(); // End time of actual VSMG function and error calculations
	auto ellapse2= (fin2-start2);
	cout <<" Run time for all Vcycle iterations = "<<ellapse2.count()/1000000000.0<<endl;

Writesolution(x);//write the solution.txt file (function is in gssolve.hpp file)

}
else if (iter==0) {
	int k=1;
	er=1.0;
	auto start2= std::chrono::system_clock::now(); //wall clock
	while (er>0.00000000001){//running Vcycle until convergence
		Vcycle(lvl,x,b,stencil);
		er =Residualnorm(x,b,stencil);
		k++;
	}
	auto fin2= std::chrono::system_clock::now(); // End time of actual VSMG function and error calculations
	auto ellapse2= (fin2-start2);
	cout <<" Run time for "<<k<<" Vcycle iterations = "<<ellapse2.count()/1000000000.0<<endl;
	cout<<"Residual norm = "<<er<<endl;

// To facilitate testing, outputs test results to file.
    ofstream reportfile;
    reportfile.open ("report.txt",std::ios_base::app);
    //errorfile << setw(20)<< "run time" << setw(20)<< "Residual norm"<< setw(20)<<"Number of levels"<<setw(20)<<"Error norm"<<endl;
    reportfile <<setw(20)<<ellapse2.count()/1000000000.0<<setw(20)<<er<< setw(20)<<lvl+1<<setw(20)<< Errornorm(x) <<endl;
    reportfile.close();
// can be commented out for testing

Writesolution(x);//write the solution.txt file (function is in gssolve.hpp file)
 

}

return 0;
}
