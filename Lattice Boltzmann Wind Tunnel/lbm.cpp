#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <map>

#include "parsing.hpp"

enum Direction {D_C,D_N,D_E,D_S,D_W,D_NE,D_SE,D_SW,D_NW};


using namespace std;

// Creates iterators for each flag
std::vector<std::vector<int>> Domain(const Grid& cells, int flag) {
        std::vector<std::vector<int>> filtered_cells;
        // Store only the cells that match the flag
        for (int j =0;j<cells.nyp();j++){
                for (int i = 0; i<cells.nxp();i++)
                 {
                        if (cells(i,j).flag == flag) {
                            filtered_cells.push_back({i,j});
                            }
                 }
        }
return filtered_cells;
};

//Pulls information from old lattice to new one, for fluid domain (given the set of indexes of all fluid cells)
void Pull(Grid &lattice,const Grid & lattice_old,const std::vector<std::vector<int>>& domain){
    #pragma omp parallel for
    for (auto index : domain ){
        lattice(index[0],index[1]).boltz[D_N] = lattice_old(index[0],index[1]-1).boltz[D_N];
        lattice(index[0],index[1]).boltz[D_E] = lattice_old(index[0]-1,index[1]).boltz[D_E];
        lattice(index[0],index[1]).boltz[D_S] = lattice_old(index[0],index[1]+1).boltz[D_S];
        lattice(index[0],index[1]).boltz[D_W] = lattice_old(index[0]+1,index[1]).boltz[D_W];
        lattice(index[0],index[1]).boltz[D_NE] = lattice_old(index[0]-1,index[1]-1).boltz[D_NE];
        lattice(index[0],index[1]).boltz[D_SE] = lattice_old(index[0]-1,index[1]+1).boltz[D_SE];
        lattice(index[0],index[1]).boltz[D_SW] = lattice_old(index[0]+1,index[1]+1).boltz[D_SW];
        lattice(index[0],index[1]).boltz[D_NW] = lattice_old(index[0]+1,index[1]-1).boltz[D_NW];
    }
}

// Pulls information into the walls and calculates the drag and lift coefficients(given the set of indexes of all bounce back cells)
std::vector<double> PullNoSlip(Grid & lattice,const std::vector<std::vector<int>> &domain){
        std::vector<double> Forces(2,0);
        double fx=0,fy=0;
        #pragma omp parallel for reduction(+:fx,fy)
        for (auto index : domain ){
                if (index[1]==0){//lower wall
                        if (index[0]==0){//left corner
                    lattice(index[0],index[1]).boltz[D_N] = lattice(index[0],index[1]+1).boltz[D_S];
                    lattice(index[0],index[1]).boltz[D_NE] = lattice(index[0]+1,index[1]+1).boltz[D_SW];
                }
                 else if (index[0]==lattice.nxp()-1){//right cornner
                    lattice(index[0],index[1]).boltz[D_N] = lattice(index[0],index[1]+1).boltz[D_S];
                    lattice(index[0],index[1]).boltz[D_NW] = lattice(index[0]-1,index[1]+1).boltz[D_SE];

                 }
                 else{
                    lattice(index[0],index[1]).boltz[D_N] = lattice(index[0],index[1]+1).boltz[D_S];
                    lattice(index[0],index[1]).boltz[D_NE] = lattice(index[0]+1,index[1]+1).boltz[D_SW];
                    lattice(index[0],index[1]).boltz[D_NW] = lattice(index[0]-1,index[1]+1).boltz[D_SE];
                        }
                 }
                else if (index[1]==lattice.nyp()-1){//top wall
                    if (index[0]==0){ //Left Corner
                    lattice(index[0],index[1]).boltz[D_S] = lattice(index[0],index[1]-1).boltz[D_N];
                    lattice(index[0],index[1]).boltz[D_SE] = lattice(index[0]+1,index[1]-1).boltz[D_NW];
                }
                else if (index[0]==lattice.nxp()-1){//right cornner
                    lattice(index[0],index[1]).boltz[D_S] = lattice(index[0],index[1]-1).boltz[D_N];
                    lattice(index[0],index[1]).boltz[D_SW] = lattice(index[0]-1,index[1]-1).boltz[D_NE];
                }
                else {
                    lattice(index[0],index[1]).boltz[D_S] = lattice(index[0],index[1]-1).boltz[D_N];
                    lattice(index[0],index[1]).boltz[D_SW] = lattice(index[0]-1,index[1]-1).boltz[D_NE];
                    lattice(index[0],index[1]).boltz[D_SE] = lattice(index[0]+1,index[1]-1).boltz[D_NW];
            }
                }
                 else { //for obstacle
                    //S
                    if (lattice(index[0],index[1]-1).flag==0) {
                            fy+=2* lattice(index[0],index[1]-1).boltz[D_N];
                            lattice(index[0],index[1]).boltz[D_S] = lattice(index[0],index[1]-1).boltz[D_N];
                    }
                    //N
                    if (lattice(index[0],index[1]+1).flag==0) {
                            fy-=2* lattice(index[0],index[1]+1).boltz[D_S];
                            lattice(index[0],index[1]).boltz[D_N] = lattice(index[0],index[1]+1).boltz[D_S];
                    }
                    //W
                    if (lattice(index[0]-1,index[1]).flag==0) {
                            fx+=2* lattice(index[0]-1,index[1]).boltz[D_E];
                            lattice(index[0],index[1]).boltz[D_W] = lattice(index[0]-1,index[1]).boltz[D_E];
                    }
                    //E
                    if (lattice(index[0]+1,index[1]).flag==0){
                            fx-=2* lattice(index[0]+1,index[1]).boltz[D_W];
                            lattice(index[0],index[1]).boltz[D_E] = lattice(index[0]+1,index[1]).boltz[D_W];
                    }
                    //SW
                    if (lattice(index[0]-1,index[1]-1).flag==0) {
                            fx+=2* lattice(index[0]-1,index[1]-1).boltz[D_NE];
                            fy+=2* lattice(index[0]-1,index[1]-1).boltz[D_NE];
                            lattice(index[0],index[1]).boltz[D_SW] = lattice(index[0]-1,index[1]-1).boltz[D_NE];
                    }
                    //SE
                    if (lattice(index[0]+1,index[1]-1).flag==0) {
                            fx-=2* lattice(index[0]+1,index[1]-1).boltz[D_NW];
                            fy+=2* lattice(index[0]+1,index[1]-1).boltz[D_NW];
                            lattice(index[0],index[1]).boltz[D_SE] = lattice(index[0]+1,index[1]-1).boltz[D_NW];
                    }
                    //NW
                    if (lattice(index[0]-1,index[1]+1).flag==0) {
                            fx+=2* lattice(index[0]-1,index[1]+1).boltz[D_SE];
                            fy-=2* lattice(index[0]-1,index[1]+1).boltz[D_SE];
                            lattice(index[0],index[1]).boltz[D_NW] = lattice(index[0]-1,index[1]+1).boltz[D_SE];
                    }
                    //NE
                    if (lattice(index[0]+1,index[1]+1).flag==0) {
                            fx-=2* lattice(index[0]+1,index[1]+1).boltz[D_SW];
                            fy-=2* lattice(index[0]+1,index[1]+1).boltz[D_SW];
                            lattice(index[0],index[1]).boltz[D_NE] = lattice(index[0]+1,index[1]+1).boltz[D_SW];
                    }
            }
        }
        Forces={fx,fy};
        return Forces;
    }

//velocity calculator(given the set of indexes of all fluid cells)
void VelocityCalculator(Grid &lattice,const std::vector<std::vector<int>>& domain){
    #pragma omp parallel for
    for (auto index : domain) {
            lattice(index[0],index[1]).velocity[0]= lattice(index[0],index[1]).boltz[D_E]+lattice(index[0],index[1]).boltz[D_SE]+lattice(index[0],index[1]).boltz[D_NE]
                                                -lattice(index[0],index[1]).boltz[D_W]-lattice(index[0],index[1]).boltz[D_SW]-lattice(index[0],index[1]).boltz[D_NW];
            lattice(index[0],index[1]).velocity[1]= lattice(index[0],index[1]).boltz[D_N]+lattice(index[0],index[1]).boltz[D_NE]+lattice(index[0],index[1]).boltz[D_NW]
                                                -lattice(index[0],index[1]).boltz[D_S]-lattice(index[0],index[1]).boltz[D_SE]-lattice(index[0],index[1]).boltz[D_SW];
                                                }
}

//Calculates the inlet boundary condition always in East direction
void InletInitializer(Grid &lattice, std::vector<std::vector<int>> &domain, double uin){
double ue=(6.0/9)*uin;
double ude=(6.0/36)*uin;
#pragma omp parallel for
for (auto index : domain){
                lattice(index[0],index[1]).boltz[D_E] = lattice(index[0]+1,index[1]).boltz[D_W]+ue;
                lattice(index[0],index[1]).boltz[D_SE] = lattice(index[0]+1,index[1]-1).boltz[D_NW]+ude;
                lattice(index[0],index[1]).boltz[D_NE] = lattice(index[0]+1,index[1]+1).boltz[D_SW]+ude;
}
}

//Calculates the outlet boundary condition
void OutletInitializer(Grid &lattice, std::vector<std::vector<int>> &domain){

#pragma omp parallel for
for (auto index : domain){
    double u2C=lattice(index[0]-1,index[1]).velocity[0]*lattice(index[0]-1,index[1]).velocity[0]+
                lattice(index[0]-1,index[1]).velocity[1]*lattice(index[0]-1,index[1]).velocity[1];

    double u2N=lattice(index[0]-1,index[1]-1).velocity[0]*lattice(index[0]-1,index[1]-1).velocity[0]+
                lattice(index[0]-1,index[1]-1).velocity[1]*lattice(index[0]-1,index[1]-1).velocity[1];

    double u2S=lattice(index[0]-1,index[1]+1).velocity[0]*lattice(index[0]-1,index[1]+1).velocity[0]+
                lattice(index[0]-1,index[1]+1).velocity[1]*lattice(index[0]-1,index[1]+1).velocity[1];

    double ux=lattice(index[0]-1,index[1]).velocity[0];
    double ucNE=lattice(index[0]-1,index[1]-1).velocity[0]+lattice(index[0]-1,index[1]-1).velocity[1];
    double ucSE=lattice(index[0]-1,index[1]+1).velocity[0]-lattice(index[0]-1,index[1]+1).velocity[1];

        lattice(index[0],index[1]).boltz[D_W] = -lattice(index[0]-1,index[1]).boltz[D_E]+2*(1.0/9)*(1+(9.0/2)*ux*ux-(3.0/2)*u2C);
        lattice(index[0],index[1]).boltz[D_SW] = -lattice(index[0]-1,index[1]-1).boltz[D_NE]+2*(1.0/36)*(1+(9.0/2)*ucNE*ucNE-(3.0/2)*u2N);
        lattice(index[0],index[1]).boltz[D_NW] = -lattice(index[0]-1,index[1]+1).boltz[D_SE]+2*(1.0/36)*(1+(9.0/2)*ucSE*ucSE-(3.0/2)*u2S);

}
}

//Updates the boltzman values, based on omega which is 1/relaxation time
void Collision(Grid &lattice,const std::vector<std::vector<int>> &domain,double omega){
 #pragma omp parallel for
 for (auto index : domain){
    std::vector<double> feq(9);
    double u2=lattice(index[0],index[1]).velocity[0]*lattice(index[0],index[1]).velocity[0]+
                lattice(index[0],index[1]).velocity[1]*lattice(index[0],index[1]).velocity[1];
    double ux=lattice(index[0],index[1]).velocity[0];
    double uy=lattice(index[0],index[1]).velocity[1];
    double ucNE=lattice(index[0],index[1]).velocity[0]+lattice(index[0],index[1]).velocity[1];
    double ucNW=-lattice(index[0],index[1]).velocity[0]+lattice(index[0],index[1]).velocity[1];
    double ucSE=lattice(index[0],index[1]).velocity[0]-lattice(index[0],index[1]).velocity[1];
    double ucSW= -lattice(index[0],index[1]).velocity[0]-lattice(index[0],index[1]).velocity[1];
    double density =lattice.Density(index[0],index[1]);
    double scale1=1.0/9;
    double scale2 =1.0/36;
    double scale3 =9.0/2;
    double scale4=3.0/2;
    //C
    feq[0]=(4.0/9) *  (density-scale4*u2);
    //N
    feq[1]=scale1 * (density+ 3*uy +scale3*(uy*uy)-scale4*u2 );
    //E
    feq[2]=scale1 *  (density+ 3*ux +scale3*(ux*ux)-scale4*u2 );
    //S
    feq[3]=scale1 * (density- 3*uy + scale3*(uy*uy)-scale4*u2 );
    //W
    feq[4]=scale1*  (density- 3*ux + scale3*(ux*ux)-scale4*u2 );
    //NE
    feq[5]=scale2 * ((density+ 3*ucNE + scale3*(ucNE*ucNE))-scale4*u2 );
    //SE
    feq[6]=scale2 * ((density+ 3*ucSE + scale3*(ucSE*ucSE))-scale4*u2 );
    //SW
    feq[7]=scale2* ((density+3*ucSW + scale3*(ucSW*ucSW))-scale4*u2 );
    //NW
    feq[8]=scale2* ((density+ 3*ucNW + scale3*(ucNW*ucNW))-scale4*u2 );

    for (int i=0;i<9;i++){
        lattice(index[0],index[1]).boltz[i]-=(omega)* (lattice(index[0],index[1]).boltz[i]-feq[i]);
    }
 }
 }


//Main program
int main(int argc, char**argv){

Grid lattice(1,1);
int timesteps,Vtksteps,DLPrint_interval=100;
std::string filename;
double omega, uin=101;
if (argc!=2){
        cout<< "Please provide the name of the parameter file for initialization\n";
        std::string paramfile;
        cin>>paramfile ;
        lattice = InputParsing(paramfile,timesteps,Vtksteps,uin,omega,filename,DLPrint_interval);

}
else { //reading param file and generating the domain accordingly

lattice = InputParsing(argv[1],timesteps,Vtksteps,uin,omega,filename,DLPrint_interval);
}
if (lattice.nxp()==1) return 1;

cout<< "Problem domain "<<lattice.nxp()-2<<" X "<<lattice.nyp()-2<< "\nTimesteps = "<<timesteps<<"\nPrint interval = " <<Vtksteps<< "\nuin = "<< uin<<endl;
/*
for (int j=0;j<lattice.nyp();j++){
    for (int i=0;i<lattice.nxp();i++){
        cout<<setw(5)<<lattice(i,j).velocity[0]<<" ";
    }
    cout<< endl;
}*/

//place holder for drag and lift respectively
std::vector<double> ForcesDL(2);

//identifying the domains for easier iteration
std::vector<std::vector<int>> fluiddomain= Domain (lattice,0);
std::vector<std::vector<int>> noslip=Domain (lattice,1);
std::vector<std::vector<int>> inlet=Domain (lattice,2);
std::vector<std::vector<int>> outlet=Domain(lattice,3);

//Boudary conditions
InletInitializer(lattice,inlet,uin);
OutletInitializer(lattice,outlet);
ForcesDL=PullNoSlip(lattice,noslip);



//second instance of the domain
Grid new_lattice=lattice;

int k=0;
int filenumber=0;
while (k < timesteps){

//calculating velocity of the fluid domain
VelocityCalculator(lattice,fluiddomain);
if (k%Vtksteps==0){
    write(filename,filenumber,lattice);
    filenumber++;
}

//collision calculations
Collision(lattice,fluiddomain,omega);

//Reinitializing boundary and walls
InletInitializer(lattice,inlet,uin);
OutletInitializer(lattice,outlet);
ForcesDL=PullNoSlip(lattice,noslip);
if (k%DLPrint_interval==0){
cout<< "\nT = "<<k <<endl;
cout<< "Drag = " <<ForcesDL[0] <<"\nLift = " << ForcesDL[1]<<endl;

}

//Pulling to new lattice
Pull(new_lattice,lattice,fluiddomain);

Dataswap(new_lattice,lattice);
//Divergence check
if (lattice(1,1).boltz[3]<0 ){
    cout << "Divergence detected" <<endl;
    break;
}
k++;
}

//Final iteration velocity calculation
VelocityCalculator(lattice,fluiddomain);

write(filename,k,lattice);
cout<< "Final iteration "<<k<<" written to file.\n";



return 0;
}
