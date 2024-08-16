#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>

using namespace std;

//multiplies a sparse matrix (in map format) to a vector
vector<double> matrix_vector_multiplication(const vector<map<int,double>> &A,const vector<double> &x) {
vector <double> result (x.size());

#pragma omp parallel for
for (long unsigned int i=0;i<A.size();i++){
    double sum=0;
    for(const auto &j : A[i]){
        sum+= x[j.first] * j.second;
    }
    result[i]=sum;
}
return result;
}

//residual norm of a system of equations with a map coefficient matrix
double RNmap(const vector<map<int,double>> &A, const vector<double> &u, const vector<double> &b){
double res=0;
#pragma omp parallel for reduction(+:res)
for (long unsigned int i=0;i<A.size();i++){
    double sum=0;
    for(const auto &j : A[i]){
        sum+= u[j.first] * j.second; //j.first is the index of the vertex and j.second is its value
    }
    sum=sum-b[i];
    res+=sum*sum;
}
return sqrt(res/A.size());
}

// uses euclidean norm to normalize the given vector
void euclidean_normalize(vector<double> &x){
double sum=0;
#pragma omp parallel for reduction (+:sum)
for (long unsigned int i=0;i<x.size();i++){
    sum+=x[i]*x[i];
}
sum = sqrt(sum);
#pragma omp parallel for
for (long unsigned int i=0;i<x.size();i++){
    x[i]=x[i]/sum;
}
}

//Gauss Seidel on a map coefficient matrix, until residual norm reaches eppsilon
void GSmap(vector<map<int,double>> &A,vector<double> &u,const vector<double> &f, double epsilon) {

double residualnorm=1;
while (residualnorm>epsilon){
//red vertices
for (long unsigned int i=0;i<A.size();i++){
    double sum=0;
    for(const auto &j : A[i]){
        sum+= u[j.first] * j.second;
    }
    sum=sum-u[i]*A[i][i];
    u[i]=(f[i]-sum)/A[i][i];
}
residualnorm = RNmap(A,u,f);//calculate the new residual norm
}
}

//dot product of 2 vectors
double vector_vector_mult(const vector<double>&a,const vector<double>&b){
double x=0;
#pragma omp parallel for reduction (+:x)
for (long unsigned int i=0;i<b.size();i++){
    x+=a[i]*b[i];
}
return x;
}

//Inverse power algorithm for map stiffness and mass matrices
//Uses GS with epsilon residual norm as the solver
double Inversepower (vector<map<int,double>>&A,vector<map<int,double>> &M,vector<double> &uh, double epsilon)  {
int k=0;
double lamda,lamdaold;
//arbitrary initial estimation of lamda
lamdaold=1;
lamda=0.5;

while (fabs((lamda-lamdaold)/lamdaold)>1/10000000000){ //fabs for absolute value of a double precision number
lamdaold=lamda;
k++;
vector<double> f=matrix_vector_multiplication(M,uh);
GSmap(A,uh,f,epsilon); //running GS solver
euclidean_normalize(uh); //normalizing uh
lamda = (vector_vector_mult(uh,matrix_vector_multiplication(A,uh)))/(vector_vector_mult(uh,matrix_vector_multiplication(M,uh)));
cout << "Iteration " << k<< " done, new estimation of lamda is " <<setprecision (15) << lamda<<endl;
}
return lamda;
}
