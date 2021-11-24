/*
First simple test
1D homogeneous problem Au = 0
Finest grid n = 64 intervals (65 points/nodes)
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
    int n = 64;
    int l1 = n+1; // fine grid
    int l2 = (l1+1)/2; // coarse grid
    double epsilon = 1e-5; // convergence criteria

    // define the coefficient matrices
    double **A1 = three_stencil(l1); // fine grid
    double **A2 = three_stencil(l2); // coarse grid

    // define the right hand side vector f1 for the fine grid
    double *f1 = new double[l1];
    initialize_vec(f1, l1);

    // define the initial guess v1 for the fine grid
    double *v1 = new double[l1];
    for(int i = 0; i < l1; i++){
        v1[i] = 0.5*(sin(3*M_PI*i/(l1-1))+sin(10*M_PI*i/(l1-1)));
    }

    // allocate residual
    double *r1 = new double[l1]; // fine grid
    double *r2 = new double[l2]; // coarse grid

    // allocate and initialize error vector
    double *e1 = new double[l1]; // fine grid
    initialize_vec(e1, l1);
    double *e2 = new double[l2]; // coarse grid
    initialize_vec(e2, l2);

    // compute the residual r1 for the fine grid
    r1 = getResidual(A1, f1, v1, l1);

    // print maximum norm of the residual
    double r_max = norm_max(r1, l1);
    cout<<"Maximum norm of residual: ";
    cout<<r_max<<endl;

    int num_iter = 1; // number of V-cycle iterations

    while(r_max > epsilon){
    //for(int i = 0; i < 2; i++){
        // first relaxation/smoothing
        v1 = GS(l1, 1, A1, f1, v1);
        r1 = getResidual(A1, f1, v1, l1);

        // restriction
        // compute the residual r1 for the coarse grid
        r2 = restrict(r1, l1);

        // second relaxation/smoothing
        initialize_vec(e2, l2);

        // solve the error equation Ae=r using GS solver 
        // or relaxation handreds of times
        e2 = GS(l2, 1, A2, r2, e2);
        double *temp = getResidual(A2, r2, e2, l2);
        double error = norm_max(temp, l2);
        double conv = 1e-7;
        while (error > conv){
            e2 = GS(l2, 1, A2, r2, e2);
            temp = getResidual(A2, r2, e2, l2);
            error = norm_max(temp, l2);
        }

        //cout<<"f1:"<<endl;
        //print_v(f1, l1);

        // prolongation
        e1 = prolong(e2, l2);

        // correct the solution v1 for the fine grid
        v1 = add_vv(v1, e1, l1);

        // third relaxation/smoothing
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
