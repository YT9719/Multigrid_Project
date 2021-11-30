/*
Recursive function
2D homogeneous problem Au = 0
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

int main(){
    int l = 4; // 
    int row = l + 1; // number of rows
    int col = l + 1; // number of columns
    int num = row * col; // number of nodes
    double epsilon = 1e-5; // convergence criteria
    int num_level = 3; // number of grids
    int pre  = 1; // number of pre-relaxation
    int post = 1; // number of post-relaxation

    // --------------allocate meomories----------------//
    double **A; // the coefficient matrix for the fine grid
    double *f  = new double[num]; // the right hand side vector
    double *v  = new double[num]; // the aproximated solution
    double *r  = new double[num]; // the residual for the fine grid

    // initialize the coefficient matrix
    A = five_stencil(row, col);
    // initialize the right hand side
    initialize_vec(f, num);
    // initialize the initial guess
    for(int i = 0; i < num; i++){
        v[i] = 0.5*(sin(3*M_PI*i/(num-1))+sin(10*M_PI*i/(num-1)));
    }

    v = GS(num, 10, A, f, v);

    int row_c, col_c, num_c;
    row_c = (row + 1) / 2;
    col_c = (col + 1) / 2;
    num_c = row_c * col_c;

    // restriction matrix
    double **R = allocate_mat(num_c, num);
    initialize_mat(R, num_c, num);
    int m, n;
    for(int i = 0; i < row_c; i++){
        for(int j = 0; j < col_c; j++){
            m = i * col_c + j;
            n = i * 2 * col + j * 2;
            if(i == 0){
                if (j == 0){
                    R[m][n]       = 1.0/4;
                    R[m][n+1]     = 1.0/8;
                    R[m][n+col]   = 1.0/8;
                    R[m][n+col+1] = 1.0/16;
                } else if(j == col_c - 1){
                    R[m][n-1]     = 1.0/8;
                    R[m][n]       = 1.0/4;
                    R[m][n+col-1] = 1.0/16;
                    R[m][n+col]   = 1.0/8;
                }else{
                    R[m][n-1]     = 1.0/8;
                    R[m][n]       = 1.0/4;
                    R[m][n+1]     = 1.0/8;
                    R[m][n+col-1] = 1.0/16;
                    R[m][n+col]   = 1.0/8;
                    R[m][n+col+1] = 1.0/16;
                }
            } else if (i == row_c - 1){
                if (j == 0){
                    R[m][n-col]   = 1.0/8;
                    R[m][n-col+1] = 1.0/16;
                    R[m][n]       = 1.0/4;
                    R[m][n+1]     = 1.0/8;
                } else if(j == col_c - 1){
                    R[m][n-col-1] = 1.0/16;
                    R[m][n-col]   = 1.0/8;
                    R[m][n-1]     = 1.0/8;
                    R[m][n]       = 1.0/4;
                }else{
                    R[m][n-col-1] = 1.0/16;
                    R[m][n-col]   = 1.0/8;
                    R[m][n-col+1] = 1.0/16;
                    R[m][n-1]     = 1.0/8;
                    R[m][n]       = 1.0/4;
                    R[m][n+1]     = 1.0/8;
                }
            }else{
                if (j == 0){
                    R[m][n-col]   = 1.0/8;
                    R[m][n-col+1] = 1.0/16;
                    R[m][n]       = 1.0/4;
                    R[m][n+1]     = 1.0/8;
                    R[m][n+col]   = 1.0/8;
                    R[m][n+col+1] = 1.0/16;
                } else if(j == col_c - 1){
                    R[m][n-col-1] = 1.0/16;
                    R[m][n-col]   = 1.0/8;
                    R[m][n-1]     = 1.0/8;
                    R[m][n]       = 1.0/4;
                    R[m][n+col-1] = 1.0/16;
                    R[m][n+col]   = 1.0/8;
                }else{
                    R[m][n-col-1] = 1.0/16;
                    R[m][n-col]   = 1.0/8;
                    R[m][n-col+1] = 1.0/16;
                    R[m][n-1]     = 1.0/8;
                    R[m][n]       = 1.0/4;
                    R[m][n+1]     = 1.0/8;
                    R[m][n+col-1] = 1.0/16;
                    R[m][n+col]   = 1.0/8;
                    R[m][n+col+1] = 1.0/16;
                }
            }
        }
    }

    // prolongation matrix
    double **P = allocate_mat(num, num_c);
    initialize_mat(P, num, num_c);
    //int m, n;
    for(int i = 0; i < row_c; i++){
        for(int j = 0; j < col_c; j++){
            m = i * col_c + j;
            n = i * 2 * col + j * 2;
            if(i == 0){
                if (j == 0){
                    P[n][m]       = 1.0/1;
                    P[n+1][m]     = 1.0/2;
                    P[n+col][m]   = 1.0/2;
                    P[n+col+1][m] = 1.0/4;
                } else if(j == col_c - 1){
                    P[n-1][m]     = 1.0/2;
                    P[n][m]       = 1.0/1;
                    P[n+col-1][m] = 1.0/4;
                    P[n+col][m]   = 1.0/2;
                }else{
                    P[n-1][m]     = 1.0/2;
                    P[n][m]       = 1.0/1;
                    P[n+1][m]     = 1.0/2;
                    P[n+col-1][m] = 1.0/4;
                    P[n+col][m]   = 1.0/2;
                    P[n+col+1][m] = 1.0/4;
                }
            } else if (i == row_c - 1){
                if (j == 0){
                    P[n-col][m]   = 1.0/2;
                    P[n-col+1][m] = 1.0/4;
                    P[n][m]       = 1.0/1;
                    P[n+1][m]     = 1.0/2;
                } else if(j == col_c - 1){
                    P[n-col-1][m] = 1.0/4;
                    P[n-col][m]   = 1.0/2;
                    P[n-1][m]     = 1.0/2;
                    P[n][m]       = 1.0/1;
                }else{
                    P[n-col-1][m] = 1.0/4;
                    P[n-col][m]   = 1.0/2;
                    P[n-col+1][m] = 1.0/4;
                    P[n-1][m]     = 1.0/2;
                    P[n][m]       = 1.0/1;
                    P[n+1][m]     = 1.0/2;
                }
            }else{
                if (j == 0){
                    P[n-col][m]   = 1.0/2;
                    P[n-col+1][m] = 1.0/4;
                    P[n][m]       = 1.0/1;
                    P[n+1][m]     = 1.0/2;
                    P[n+col][m]   = 1.0/2;
                    P[n+col+1][m] = 1.0/4;
                } else if(j == col_c - 1){
                    P[n-col-1][m] = 1.0/4;
                    P[n-col][m]   = 1.0/2;
                    P[n-1][m]     = 1.0/2;
                    P[n][m]       = 1.0/1;
                    P[n+col-1][m] = 1.0/4;
                    P[n+col][m]   = 1.0/2;
                }else{
                    P[n-col-1][m] = 1.0/4;
                    P[n-col][m]   = 1.0/2;
                    P[n-col+1][m] = 1.0/4;
                    P[n-1][m]     = 1.0/2;
                    P[n][m]       = 1.0/1;
                    P[n+1][m]     = 1.0/2;
                    P[n+col-1][m] = 1.0/4;
                    P[n+col][m]   = 1.0/2;
                    P[n+col+1][m] = 1.0/4;
                }
            }
        }
    }

    print_m(P, num, num_c);

}