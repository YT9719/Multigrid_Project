/*
First simple test
1D homogeneous problem Au = 0
Finest grid n = 64 points
Two level grid 
The initial guess is two sine waves (low frequency error + high frequency error)
v_0(x) = 0.5*(sin(3*M_PI*x/n)+sin(10*M_PI*x/n))
*/

#define _USE_MATH_DEFINES //to include the mathematical constants
#include <iostream>
#include <iomanip> // include the setprecision funtion
#include <cmath> // include math functions
#include "../multigrid.h" // include all modules and functions

using namespace std;

int main(){
    int n = 8;
    int l1 = n; // fine grid
    int l2 = (n+1)/2; // coarse grid

    // define the coefficient matrix A
    double **A = allocate_mat(l1, l1);
    initialize_mat(A, l1, l1);
    for(int i = 0; i < l1; i++){
        if(i == 0){
            A[i][i] = 1; A[i][i+1] = 0.5;
        }
        else if(i == l1 - 1){
            A[i][i] = 1; A[i][i-1] = 0.5; 
        }
        else{
            A[i][i] = 1; A[i][i+1] = 0.5; A[i][i-1] = 0.5;
        }
    }

    // define the right hand side vector
    double *f = new double[l1];
    for(int i = 0; i < l1; i++){
        f[i] = 0;
    }

    // define the initial guess
    double *v = new double[l1];
    for(int i = 0; i < l1; i++){
        v[i] = 0.5*(sin(3*M_PI*i/(l1-1))+sin(10*M_PI*i/(l1+1)));
    }

    // compute the residual
    double *residual = new double[l1];
    residual = getResidual(A, f, v, l1);

    //print maximum norm of the residual
    cout<<norm_max(residual, l1)<<endl;

    // first relaxation/smoothing
    v = GS(l1, 1, A, f, v);
    residual = getResidual(A, f, v, l1);

    //print maximum norm of the residual
    cout<<norm_max(residual, l1)<<endl;

    return 1;
}
