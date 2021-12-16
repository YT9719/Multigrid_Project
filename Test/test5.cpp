/*
Parallel recursive function
2D homogeneous problem Au = 0
Finest grid row & col = 64 intervals (65^2 points/nodes)
Three level grid 
The initial guess is two sine waves (low frequency error + high frequency error)
v_0(x) = 0.5*(sin(3*M_PI*x/n)+sin(10*M_PI*x/n))
*/

#define _USE_MATH_DEFINES //to include the mathematical constants
#include <iomanip> // include the setprecision & setw funtion
#include <cmath> // include math functions
#include "../Src/multigrid_2D_parallel.h" // include all modules and functions

using namespace std;

// -----------------input data---------------------//
int n = 64; // number of intervals
int row = n - 1; // number of rows
int col = n - 1; // number of columns
int num = row * col; // number of nodes
double epsilon = 1e-6; // convergence criteria
int num_level = 3; // number of grids
int pre  = 2; // number of pre-relaxation
int post = 2; // number of post-relaxation
int n_thread = 6; // number of threads
// ------------------------------------------------//

int main(){
    omp_set_num_threads(n_thread);
    int level = 1;

    // --------------allocate meomories----------------//
    double **A; // the coefficient matrix for the fine grid
    double *f  = new double[num]; // the right hand side vector
    double *v  = new double[num]; // the aproximated solution
    double *r  = new double[num]; // the residual for the fine grid

    double dx = 1.0/(row+1);

    // initialize the coefficient matrix
    A = five_stencil(row, col, dx);
    // initialize the right hand side
    initialize_vec(f, num);
    // initialize the initial guess
    for(int i = 0; i < num; i++){
        v[i] = 0.5*(sin(3*M_PI*(i+1)/(num+1))+sin(10*M_PI*(i+1)/(num+1)));
    }
    // initialize the max norm of residual
    double r_max = 1;
    // initialize the number of iterstions
    int iter = 1;

    double start_time = omp_get_wtime();
    // print the titles
    cout<<setw(10)<<"Iterations"<<"   "<<"Residual"<<endl;

    while(r_max > epsilon){

        // multigrid V-cycle
        v = V_cycle(A, v, f, row, col, level, pre, post, num_level);

        // compute the residual for the fine grid
        residual(r, A, f, v, num);

        // compute maximum norm of the residual
        r_max = norm_max(r, num);

        // print the number of itertions and max norm of the residual
        cout<<setw(10)<<iter<<" "<<setprecision(3)<<r_max<<endl;
        iter = iter + 1;
    }

    /*
    cout<<"Approximated solution:"<<endl;
    for(int i = 0; i < num; i++){
        cout<<fixed<<setprecision(5)<<v[i]<<endl;
    }
    */

    double run_time = omp_get_wtime() - start_time;
    cout<<"time:"<<run_time<<endl;

    return 1; 
}