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
#include "../Src/multigrid_1D.h" // include all modules and functions
#include <chrono>
using namespace std::chrono;

using namespace std;

// -----------------input data---------------------//
int n = 64; // number of intervals
int num = n + 1; // number of nodes
double epsilon = 1e-7; // convergence criteria
int num_level = 3; // number of grids
int pre  = 1; // number of pre-relaxation
int post = 1; // number of post-relaxation
// ------------------------------------------------//

int main(){
    auto start = high_resolution_clock::now();  
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
    // initialize the max norm of residual
    double r_max = 1;
    // initialize the number of iterstions
    int iter = 1;

    // print the titles
    cout<<setw(10)<<"Iterations"<<"   "<<"Residual"<<endl;

    while(r_max > epsilon){

        // multigrid V-cycle
        v = F_cycle(A, v, f, num, level, pre, post, num_level);

        // compute the residual for the fine grid
        r = getResidual(A, f, v, num);
    
        // compute maximum norm of the residual
        r_max = norm_max(r, num);

        // print the number of itertions and max norm of the residual
        cout<<setw(10)<<iter<<" "<<setprecision(3)<<r_max<<endl;
        iter = iter + 1;
    }

    cout<<"Approximated solution:"<<endl;
    for(int i = 0; i < num; i++){
        cout<<fixed<<setprecision(5)<<v[i]<<endl;
    }

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
  
    cout << "Time taken by function: "
         << duration.count() << " microseconds" << endl;

    return 1; 
}

