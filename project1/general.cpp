#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include "time.h"
// use namespace for output and input
using namespace std;
using namespace arma;
// object for output files
ofstream ofile;
// Functions
inline double f(double x) {return 100.0*exp(-10.0*x);} // source term f(x)
inline double u(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);} // closed-form solution

int main(int argc, char *argv[]){
  int exponent;
  string outputfilename;

  if(argc <= 1){
    cout << "Two few arguments in " << argv[0] <<
    ". Read in name of output file and the exponent of the maximum number of gridpoints" << endl;
    exit(1);
  }
  else{
    outputfilename = argv[1];
    exponent = atoi(argv[2]);
  }

  // Looping over different number of gridpoints 10^nsed
  for (int i = 1; i <= exponent; i++){
    int n = (int) pow(10.,i);
    string outputfile = outputfilename;
    string numbr = to_string(i);
    outputfile.append(numbr);

    vec v(n+2, fill::zeros);
    vec x(n+2); x(0) = 0.; x(n+1) = 1;

    double h =1.;
  }
}
