/*
Second simple test
Solve the same problem as test 1 only using GS solver
1D homogeneous problem Au = 0
Finest grid n = 64 intervals (65 points/nodes)
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
    int n = 64;
    int l1 = n+1; // fine grid
    int l2 = (l1+1)/2; // coarse grid
    double epsilon = 1e-5; // convergence criteria

    // define the coefficient matrix A1 for the fine grid
    double **A1 = three_stencil(l1);

    // define the right hand side vector f1 for the fine grid
    double *f1 = new double[l1];
    initialize_vec(f1, l1);

    // define the initial guess v1 for the fine grid
    double *v1 = new double[l1];
    for(int i = 0; i < l1; i++){
        v1[i] = 0.5*(sin(3*M_PI*i/(l1-1))+sin(10*M_PI*i/(l1-1)));
    }

    // compute the residual r1 for the fine grid
    double *r1 = new double[l1];
    r1 = getResidual(A1, f1, v1, l1);
    
    // print maximum norm of the residual
    double r_max = norm_max(r1, l1);
    cout<<"Maximum norm of residual: ";
    cout<<r_max<<endl;

    int num_iter = 1; // number of V-cycle iterations

    //while(r_max > epsilon){
    while(r_max > epsilon){

        // first relaxation/smoothing
        v1 = GS(l1, 1, A1, f1, v1);
        r1 = getResidual(A1, f1, v1, l1);

        // print maximum norm of the residual
        r_max = norm_max(r1, l1);
        cout<< "Number of V-cycles: "<< num_iter<<";   ";
        cout<< "Maximum norm of residual: ";
        cout<< r_max <<endl;
        num_iter = num_iter + 1;

    }
    
    // print the final solution
    //cout<<"Approximated solution:"<<endl;
    //print_v(v1, l1);

    return 1;
}