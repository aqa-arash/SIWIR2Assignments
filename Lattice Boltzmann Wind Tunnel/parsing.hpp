#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>



using namespace std;

//underlying structure which holds all information of each cell
struct Cell {

    int flag;

    double density;
    vector<double> velocity;
    vector<double> boltz; //C N E S W NE SE SW NW

    // Constructor to initialize vectors
    Cell()
        : flag(0),density(1),velocity(3), boltz(9) {}
};



// Grid class to control 2D cells
class Grid
{
public:
  // allocate memory in constructor
  Grid(int __nx__, int __ny__)
      : nx(__nx__)
      , ny(__ny__)
  {
      data = std::vector<Cell>(nx * ny); // Initialize vector with size
  }

  // free memory in destructor (not required for std::vector)
  ~Grid() {} // Empty destructor


  //swapping grids
  friend void Dataswap(Grid &a,Grid &b);


  // Quick assign method
    void CellVectorAllSet(vector<double> input,char info='v')
    {
        if (info=='v'){
        #pragma omp parallel for
        for (int j=0;j<ny;j++) {
            for (int i=0;i<nx;i++) {
                data[i + j * nx].velocity=input;
                }
        }
        }
        else if (info=='l'){
            #pragma omp parallel for
            for (int j=0;j<ny;j++) {
                for (int i=0;i<nx;i++) {
                    data[i + j * nx].boltz=input;
                }
        }
        }
    }

    //even quicker - unity
     void CellValueAllSet(double input, char type ='d')
    {
        if (type=='d'){
        #pragma omp parallel for
        for (int k=0;k<ny*nx;k++) {
           data[k].density=input;
        }
    }
    else if (type=='f'){
        #pragma omp parallel for
        for (int k=0;k<ny*nx;k++) {
           data[k].flag=int(input);
        }
    }

    }

    //Boundary set
    void TunnelFlagInitialize (){
        //left and right
        #pragma omp parallel for
            for (int j=0;j<ny;j++){
                    data[j*nx].flag=2;
                    data[(j)*nx+nx-1].flag=3;}
        //top and bottom
        #pragma omp parallel for
            for (int i=0;i<nx;i++){
                    data[i + (ny-1) * nx].flag=1;
                    data[i].flag=1;}
    }

    double Density(int i, int j){
        double d=0;
        for (int ii =0;ii<9;ii++){
            d+=data[i + j * nx].boltz[ii];
        }
        return d;
    }

  // getter
  auto operator()(int i, int j) const { return data[i + j * nx]; }

  // setter
  auto& operator()(int i, int j) { return data[i + j * nx]; }

  // access dimensions- Domain size is considered 1x1
  int nyp() const { return ny; }//number of points in y direction
  int nxp() const { return nx; }// number of points in x direction



private:
  // internal data
  int nx, ny;
  std::vector<Cell> data;
};


void Dataswap(Grid &a,Grid &b){
  using std::swap;
  swap(a.data,b.data);
  };

//writes vtk file
void write(std::string filename, int timestep, Grid & lattice){
    int Nx=lattice.nxp();
    int Ny=lattice.nyp();
        std::ofstream file(filename + std::to_string(timestep) + ".vtk");
        file << "# vtk DataFile Version 4.0\n";
        file << "SiWiRVisFile\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING 1 1 1\n";
        file << "POINT_DATA " << Nx * Ny << "\n";

        file << "SCALARS flags unsigned_int 1\n";
        file << "LOOKUP_TABLE default\n";
        for(int y=0; y<Ny; ++y){
            for(int x=0; x<Nx; ++x){
                file << lattice(x,y).flag << endl;
            }
        }

        file << "SCALARS density double 1\n";
        file << "LOOKUP_TABLE default\n";
        for(int y=0; y<Ny; ++y){
            for(int x=0; x<Nx; ++x){
                file << "1" << endl; // Write the density value after calculating density
            }
        }

        file << "VECTORS velocity double\n";
        for(int y=0; y<Ny; ++y){
            for(int x=0; x<Nx; ++x){
                file << lattice(x,y).velocity[0]<<" " <<lattice(x,y).velocity[1] <<" " <<lattice(x,y).velocity[2] <<  endl; // Write the velocity value after calculating velocity Ux and Uy
            }
        }
        file.close();
    };

//Flags for spherical obstacle in the flow path
void ObstacleFlag(Grid & A, int cx,int cy, double r){
    int searchdomain=int(r);
    for (int j=cy-searchdomain;j<cy+searchdomain;j++){
        for (int i=cx-searchdomain;i<cx+searchdomain;i++){
                double distance= (i+0.5-cx)*(i+0.5-cx)+(j+0.5-cy)*(j+0.5-cy);
                if (distance<r*r){
                    A(i,j).flag=1;
                }
    }
    }}

//returns the domain based on the input pgm file
Grid PGMparsing (std::string filename){

        cout<< "Parsing pgm Initiated \n";
            ifstream file(filename);
            string line,dataname;
            int k =0;
            int Nx=-1,Ny=-1;
            int i=0,j=0;
            bool ValuesStarted =false;
//going through the file, line by line
  while (std::getline(file, line) && !ValuesStarted) {
        std::istringstream iss(line);
            std::string val;
            int n=0;
            while (iss >> val) {
                    if (k==0 || (n==0 && val == "#")) break;
                    if (n==0) Nx=stoi(val);
                    else {
                        Ny=stoi(val);
                        ValuesStarted=true;
                        }
                n++;
                }
        k++;
  }

Grid A(Nx+2,Ny+2);
A.TunnelFlagInitialize();
A.CellValueAllSet(1);
A.CellVectorAllSet({0,0,0});
A.CellVectorAllSet({4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36},'l');

j=Ny;
  //going through the file, line by line
  while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string val;
        int temp;
        while (iss >> val) {
                temp=stoi(val);
                if (i%Nx!=0 || i==0){
                    if (temp<255) A(i+1,j).flag=1;
                    i++;
                }
                else {
                    i=0;
                    j--;
                    if (temp<255) A(i+1,j).flag=1;
                    i++;
                }
            }

  }
  return A;  }


//reads the param file and generates the domain, modifying required values
Grid InputParsing(std::string Parameterfile,int &timesteps,int &Vtksteps,double &uin,double &omega,std::string &outputfile,int &DLPrint_interval) {
            ifstream file(Parameterfile);
            string line,dataname;
            int Nx=-1,Ny=-1,spherex=-1,spherey=-1,diameter=0;
            double reinolds_number=0;
            string pgmfile="NoFile";
  while (std::getline(file, line)) {
        std::istringstream iss(line);
            std::string val;
            while (iss >> val) {

                        if (val=="sizex"){
                            iss>>val;
                            Nx=stoi(val);

                        }
                        else if (val=="sizey"){
                            iss>>val;
                            Ny=stoi(val);


                        }
                        else if (val=="timesteps"){
                            iss>>val;
                            timesteps=stoi(val);

                        }
                        else if (val=="uin"){
                            iss>>val;
                            uin=stod(val);


                        }
                        else if (val=="Re"){
                            iss>>val;
                            reinolds_number=stod(val);

                        }
                        else if (val=="spherex"){
                            iss>>val;
                            spherex=stoi(val);

                        }
                        else if (val=="spherey"){
                            iss>>val;
                            spherey=stoi(val);


                        }
                        else if (val=="diameter"){
                            iss>>val;
                            diameter=stoi(val);


                        }
                        else if (val=="vtk_file"){
                            iss>>val;
                            outputfile=val;


                        }
                        else if (val=="vtk_step"){
                            iss>>val;
                            Vtksteps=stoi(val);

                        }
                        else if (val== "geometry"){
                            iss>>val;
                            pgmfile=val;
                            }
                        else if (val== "Print_Interval"){
                            iss>>val;
                            DLPrint_interval=stoi(val);
                            }
                }
            }

if (pgmfile!="NoFile"){
Grid lattice= PGMparsing(pgmfile);
double viscosity=uin*(lattice.nyp()-2)/reinolds_number;
double taw = 3*viscosity + 0.5;
cout<< "Based on the input variables taw is = "<< taw<<endl;
if (taw<0.51) cout<< "Consider modifying input variables, Could cause unstability! \n";
omega = 1.0/taw;
return lattice;
}
else {
double viscosity=uin*Ny/reinolds_number;
double taw = 3*viscosity + 0.5;
cout<< "Based on the input variables taw is = "<< taw<<endl;
if (taw<0.51) cout<< "Consider modifying input variables, Could cause unstability! \n";
omega = 1.0/taw;
    Grid lattice(Nx+2,Ny+2);
    lattice.TunnelFlagInitialize();
    lattice.CellValueAllSet(1);
    lattice.CellVectorAllSet({0,0,0});
    lattice.CellVectorAllSet({4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36},'l');
if (spherex!=-1 &&spherey!=-1 && diameter!=0){
    double radius = diameter/2.0;
    ObstacleFlag(lattice, spherex,spherey,radius);
    //obstacles flagged
    }
return lattice;
}
}
