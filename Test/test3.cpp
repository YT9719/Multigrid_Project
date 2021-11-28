/*
Recursive function
1D homogeneous problem Au = 0
Finest grid n = 64 intervals (65 points/nodes)
Three level grid 
The initial guess is two sine waves (low frequency error + high frequency error)
v_0(x) = 0.5*(sin(3*M_PI*x/n)+sin(10*M_PI*x/n))
*/

#define _USE_MATH_DEFINES //to include the mathematical constants
#include <iostream>
#include <iomanip> // include the setprecision funtion
#include <cmath> // include math functions
#include "../multigrid.h" // include all modules and functions

using namespace std;

// -----------------input data---------------------//
int n = 64; // number of intervals
int num = n + 1; // number of nodes
double epsilon = 1e-5; // convergence criteria
int num_level = 3; // number of grids
int pre  = 1; // number of pre-relaxation
int post = 1; // number of post-relaxation
// ------------------------------------------------//



double *V_cycle(double **A, double *v, double *f, int num, int level){
    
    // pre-relaxation
    A = three_stencil(num);
    v = GS(num, pre, A, f, v);

    // compute the residual for the fine grid
    double *r = getResidual(A, f, v, num);

    // restriction 
    double *rc = restrict(r, num);
    level = level + 1;
    int num_c = (num + 1) / 2;

    // inititalize the error vector
    double *e = new double[num_c];
    initialize_vec(e, num_c);

    // inititalize the coefficient matrix
    double **A_c = three_stencil(num_c);

    // stop recursion at coarsest grid, otherwise continue recursion
    if(level == num_level){
        //solve the error equation
        e = GS(num_c, 1, A_c, rc, e);
        double *temp = getResidual(A_c, rc, e, num_c);
        double error = norm_max(temp, num_c);
        double conv = 1e-7;
        while (error > conv){
            e = GS(num_c, 1, A_c, rc, e);
            temp = getResidual(A_c, rc, e, num_c);
            error = norm_max(temp, num_c);
        }
    } else{
        e = V_cycle(A_c, e, rc, num_c, level);
    }

    // prolongation 
    //cout<<"e:"<<endl;
    //print_v(e, num_c);
    double *ef = prolong(e, num_c);

    // correct the approximated solution
    v = add_vv(v, ef, num);

    // post-relaxation
    v = GS(num, post, A, f, v);

    // compute the residual for the fine grid
    //r = getResidual(A, f, v, num);

    return v;
}


int main(){
    int level = 1;

    // --------------allocate meomories----------------//
    double **A; // the coefficient matrix for the fine grid
    double *f  = new double[num]; // the right hand side vector
    double *v  = new double[num]; // the aproximated solution
    double *r  = new double[num]; // the residual for the fine grid

    // initialize the coefficient matrix
    A = three_stencil(num);
    // initialize the right hand side
    initialize_vec(f, num);
    // initialize the initial guess
    for(int i = 0; i < num; i++){
        v[i] = 0.5*(sin(3*M_PI*i/(num-1))+sin(10*M_PI*i/(num-1)));
    }

    for(int i = 0; i < 50; i++){

        // multigrid V-cycle
        v = V_cycle(A, v, f, num, level);

        // compute the residual for the fine grid
        r = getResidual(A, f, v, num);

        // print maximum norm of the residual
        double r_max = norm_max(r, num);
        cout<<"Maximum norm of residual: ";
        cout<<r_max<<endl;
    }

    //cout<<"v:"<<endl;
    //print_v(v, num);

    return 1; 
}






