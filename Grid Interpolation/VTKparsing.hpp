#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

void read(std::vector<std::vector<double>> &vertices,std::vector<std::vector<int>> &faces,std::vector<std::vector<double>> &values,const std::string &filename) {

  ifstream file(filename);
  string line,dataname;

  //int valuesize=1;

bool POINT_DATA=false,CELLS =false,POINTS =false, valuesstarted=false; //check which data we are reading.
//auto datatype=double;

int k =0;
//going through the file, line by line
  while (std::getline(file, line)) {
        std::istringstream iss(line);
        //parsing points
        if (POINTS){
                std::vector<double> coordinates;
                double val;
                while (iss >> val ) {
                    coordinates.push_back(val);
                    }
                vertices.push_back(coordinates);
                k--;
                if (k==0){
                    POINTS=false;
                }
           }

        //parsing cells
        else if (CELLS) {
                int n=0;
                std::vector<int> triangle;
                int val;
                while (iss >> val ) {
                    if (n>0) triangle.push_back(val);
                    n++;
                    }
                faces.push_back(triangle);
                k--;
                if (k==0){
                    CELLS=false;

                }

        }


        //parsing dara
        else if (POINT_DATA){
                int n=0;
                bool SCALARS=false;
                std::vector<double> pointvalue;

                if (!valuesstarted) {
                    string val;
                    while (iss >> val ) {

                        if (val=="SCALARS" && n==0) {SCALARS=true; n++; continue;}
                        if (SCALARS && n==1) {dataname=val; n++; continue;}
                        if (SCALARS && n==2) {n++; continue;
                            //if (val=="float" || val=="double") continue;
                            //if (val=="int") { datatype = int;  continue; }
                            }
                        if (SCALARS && n==3) { n++; continue;}

                        if (val=="LOOKUP_TABLE") break;
                        try {
                                    pointvalue.push_back(stod(val));
                                    valuesstarted=true;
                                    break;
                                    }
                        catch(...) {break;}
                            }

                        }

                if (valuesstarted){
                    double val;
                    while (iss >> val ) {
                            pointvalue.push_back(val);
                            }
                    values.push_back(pointvalue);
                    k--;
                    if (k==0) {
                        valuesstarted=false;
                        POINT_DATA=false;
                    }
                }

        }


        //parsing initial lines (not values)
        else {
                    string val;
                    int n=0;
                    while (iss >> val) {
                        if (val=="#") break;

                        if (val=="POINTS" && n==0) POINTS=true;
                        if (val=="CELLS"&& n==0) CELLS=true;
                        if (val=="POINT_DATA"&& n==0) POINT_DATA=true;

                        if (!CELLS && !POINTS && !POINT_DATA) break;

                        if (POINTS && n==1) {
                                k=stoi(val);
                        }

                        if (CELLS && n==1) {
                                k=stoi(val);
                                faces.reserve(k);
                        }


                        if (POINT_DATA && n==1) {
                                k=stoi(val);
                                //cout<<"Point data collection started with k = " << k<<endl;
                        }
                        n++;
                        }
            }

        }
 }



void write(std::vector<std::vector<double>> &vertices,std::vector<std::vector<int>> &faces,std::vector<std::vector<double>> &values, const std::string &filename) {

  std::ofstream file;
  file.open (filename);
  file << "# vtk DataFile Version 2.0\nSample triangulation\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  file << "\n";

  // Write Vertices
  file << "POINTS " << vertices.size() << " float\n";
  for (const auto& vertex : vertices){ // Iterates over the vertices vector
    for (const auto& val : vertex){ // Iterates over each element in vertex vector
        file << setprecision(13)<< val << " ";
    }
    file << "\n";
  }

  file << "\n";

  // Write Cells
  int cell_size = 0;
  for(const auto& face : faces){
    cell_size += face.size() + 1;
  }

  file << "CELLS " << faces.size()<< " " << cell_size << "\n";
  for (const auto& face : faces){
    file << face.size() << " ";
    for (const auto& val : face){
        file << val << " ";
    }
    file << "\n";
  }

  file << "\n";

  // Writing Cell types
  file << "CELL_TYPES " << faces.size() << "\n";
  for (long unsigned int i = 0 ; i < faces.size(); i++){
    file << "5\n"; // All cells are trianle i.e VTK triangle
  }

  file << "\n";

  // Write point data
  file << "POINT_DATA " << vertices.size() << "\n";
  file << "SCALARS value float " << "1\n";
  file << "LOOKUP_TABLE default\n";
  for (const auto& value : values){
    for (const auto& val : value){
        file << val ;
    }
    file << "\n";
  }

  file.close();
 }

