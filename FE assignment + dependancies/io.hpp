#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

//this code reads the text files and stores the vetex and faces data in corresponding 2D vectors
//it only works with files that are formatted exactly as the example
// you have to give it the 2D vectors and the name of the file as input arguments
void input(std::vector<std::vector<double>> &vertices,std::vector<std::vector<int>> &faces,const std::string &filename) {
std::ifstream file(filename);
std::string line;

int k =0;
while (std::getline(file, line)) {
    std::istringstream iss(line);
    if (k<3){
        std::vector<double> row;
        double val;
        int n=0;
        while (iss >> val && n<4 ) {
            n++;
                row.push_back(val);
}
            if (row.size()==3){
            vertices.push_back({row[1],row[2]});
            }
            else {k++;}
        }

    else if (k>=3){
        int val;
        std::vector<int> row;
        int n=0;
        while (iss >> val && n<4 ) {
            n++;
                row.push_back(val);
}
            if (row.size()==3){
            faces.push_back(row);
            }
            else {k++;}
        }
  }
 }

