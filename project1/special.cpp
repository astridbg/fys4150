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

  // Looping over different number of gridpoints 10^n
  for (int i = 1; i <= exponent; i++){
    int n = (int) pow(10.,i);
    string outputfile = outputfilename;
    string numbr = to_string(i);
    outputfile.append(numbr);

    // Defining the solution vector v and grid point vector x
    vec v(n+2, fill::zeros);
    vec x(n+2);
    // Defining the start and end value of the grid
    x(0) = 0.; x(n+1) = 1;
    // Defining the step length h
    double h = 1./n;
    double hh = h*h;
    // Defining the matrix element vectors of the tridiagonal matrix A
    vec a(n); vec b(n); vec c(n);
    // Defining the source term vector d = Av
    vec d(n);
    // Defining the new vectors after tridiagonal matrix algorithm
    vec c_new(n); vec d_new(n);

    // Filling vectors
    for (int i = 0; i < n; i++){
      a(i) = -1.; b(i) = 2.; c(i) = -1.;
      x(i+1) = (i+1)*h;
      d(i) = hh*f(x(i+1));
    }

    // Compute elapsed time
    clock_t start, finish; // declaring start and finish time
    start = clock(); // start time of computing

    // Forward substitution
    for(int i = 0; i < n; i++){
      c_new(i) = -((double) i+1)/(i+2);
      for(int j = 0; j <= i; j++){
        d_new(i) += ((double) j+1)/(i+2)*d(j);
      }
    }

    // Backward substitution
    v(n-1) = d_new(n-1);
    for (int i = n-2; i > 0; i--){
        v(i) = d_new(i) - c_new(i)*v(i+1);
    }

    finish = clock(); // finish time of computing
    double time = ((double)(finish - start)/CLOCKS_PER_SEC);
    cout << setprecision(32) << "Time used for " << n << " gridpoints: "<< time << endl;

    // Opening file and writing out results
    ofile.open(outputfile);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "     x:             approx:        exact:         relative error:" << endl;
    for (int i = 0; i < n+2; i++){
      double RelErr = fabs((v(i)-u(x(i)))/u(x(i))); // Compute relative error
      ofile << setw(15) << setprecision(8) << x(i);
	    ofile << setw(15) << setprecision(8) << v(i);
	    ofile << setw(15) << setprecision(8) << u(x(i));
      ofile << setw(15) << setprecision(8) << RelErr << endl;
    }
    ofile.close();
  }
  return 0;
}
