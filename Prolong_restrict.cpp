/*
Basically, the prolongation(also called interpolation) 
transfers the data on a coarse grid to a fine grid. 
The simple way is to do a linear interplation to obtain
the data on the fine grid. 

The restriction refers to transfer the data on a find 
grid to a coarse grid. The simple way is to give the 
values on the fine grid to the coarse grid. 

Right now, the code focuses on a 1D problem.
*/
#define _USE_MATH_DEFINES //to include the mathematical constants
#include <iostream>
#include <iomanip> // include the setprecision funtion
#include <cmath> // include math functions
#include "multigrid.h" // include all modules and functions

using namespace std;

int main() {
    int n=4;
    double *x1 = new double[n]; //coarse grid
    double *x2 = new double[2*n-1]; // fine grid
    double *x3 = new double[n]; //back to the coarse grid

    // assign values for the coarse grid
    for(int i=0; i<n; i++){
        x1[i] = sin(i*M_PI/n);
    }

    // print x1 vector
    cout<<"x1:"<<endl;
    print_v(x1,n);

    // apply prolongation matrix to the coarse grid
    x2 = prolong(x1, 2*n-1);

    // print x2 vector
    cout<<"x2:"<<endl;
    print_v(x2,2*n-1);

    x3 = restrict(x2, 2*n-1);
    // print x2 vector
    cout<<"x3:"<<endl;
    print_v(x3,n);

    return 1;
}

// Next step: apply prolongation and restriction to 2D problem

