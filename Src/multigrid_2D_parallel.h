// header files for 2D problem 

#include <iostream>
#include <algorithm> // include max_element function
#include <omp.h> // include OpenMP

using namespace std;

// allocate a matrix 
// row and col are number of rows and columns
double **allocate_mat(int row, int col){
    double **m = new double *[row];
    for(int i=0; i<row; i++){
        m[i] = new double [col];
    }
    return m;
}

// initialize a matrix to all zero
// row and col are number of rows and columns
void initialize_mat(double **m, int row, int col)
{
    for(int i = 0; i < row; i++) {
        for(int j = 0; j < col; j++) {
            m[i][j] = 0;
        }
    }
}

// generate coefficient matrix for five-point stencil
// row and col are number of rows and columns of the grid
double **five_stencil(int row, int col){
  int i, j, d;
  d = row * col;
  double **m = allocate_mat(d, d);
  initialize_mat(m, d, d);
  for(int k = 0; k < d; k++){
      i = k/row; j = k-i*row;
      if(i>0)     {m[k][k-row] = -0.25;}
      if(i<col-1) {m[k][k+row] = -0.25;}
      if(j>0)     {m[k][k-1] = -0.25;}
      if(j<row-1) {m[k][k+1] = -0.25;}
      m[k][k] = 1.0;
    }
  return m;
}

// initialize a vector to all zero
// n is the size of the vector
void initialize_vec(double *v, int n){
  for(int i = 0; i < n; i++){
        v[i] = 0;
    }
}

// matrix and vector multiplication
// **m is a matrix
// row and col are number of rows and columns
// *v is a vector
double *multiply_mv(double **m, int row, int col, double *v){
    double *b = new double [row];
    int i, j;
    #pragma omp parallel default(none)\
    shared(row, col, m, v, b) private(i, j)
    // private variables should not be initialized
    {   
        #pragma omp for
        for(i = 0; i < row; i++){
            double sum = 0;
            for(j = 0; j < col; j++){
                sum = sum + m[i][j] * v[j];
            }
            b[i] = sum;
        }
    }
    return b;
}

// vector addition
// *v1 & *v2 are two vectors
// n is the size of the vector
double *add_vv(double *v1, double *v2, int n){
  double *v3 = new double[n];
  for(int i = 0; i < n; i++){
    v3[i] = v1[i] + v2[i];
  }
  return v3;
}

// maximum norm of a vector
// *v is a vector
// n is the size of the vector
double norm_max(double *v, int n){
  double *v_new = new double[n]; 
  for(int i = 0; i < n; i++){
    v_new[i] = abs(v[i]);
  }
  double max = *max_element(v_new, v_new+n);
  delete[] v_new;
  return max;
}

// Restriction module（2D）
// Full-weighting restriction
// *rf is the residual on the fine grid
// row and col are number of rows and columns of fine grid
// *rc is the residual on the coarse grid
double *restrict_2D(double *rf, int row, int col){
  int num, row_c, col_c, num_c;
  num = row * col; // number of nodes of fine grid 
  row_c = (row + 1) / 2; // number of rows of coarse grid
  col_c = (col + 1) / 2; // number of columns of coarse grid
  num_c = row_c * col_c; // number of nodes of coarse grid
  double *rc = new double[num_c];
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
    rc = multiply_mv(R, num_c, num, rf);
    delete[] R;
    return rc;
}

// Prolongation module（2D）
// Full-weighting prolongation
// *ec is the error on the coarse grid
// row_c and col_c are number of rows and columns of coarse grid
// *ef is the error on the fine grid
double *prolong_2D(double *ec, int row_c, int col_c){
    int num_c, num, row, col;
    num_c = row_c * col_c; // number of nodes of coarse grid
    row = row_c * 2 - 1; // number of rows of fine grid
    col = col_c * 2 - 1; // number of columns of fine grid
    num = row * col; // number of nodes of fine grid 
    double *ef = new double [num];
    double **P = allocate_mat(num, num_c);
    initialize_mat(P, num, num_c);
    int m, n;
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
    ef = multiply_mv(P, num, num_c, ec);
    delete[] P;
    return ef;
}

// parallel Jacobi solver
// n is the szie of vector
// **M is the coefficient matrix
// *v is the right hand side vector
// *x is the initial guess of solution 
double* Jac_p(int n, int iter, double **M, double *v, double *x) {
  double sigma = 0;
  double *b = new double[n];
  for (int k = 0; k < iter; k ++) {
    int i, j;
    #pragma omp parallel default(none)\
    private(i, j, sigma) shared(n, M, x, v, b)
    {
        #pragma omp for 
        for (i = 0; i < n; i ++) {
            sigma = 0;
            for (j = 0; j < n; j ++) {
                if (j != i) {
                sigma += M[i][j]*x[j];
                }
            }
        b[i] = (v[i]-sigma)/M[i][i];
        }
    }
  }
  return b;
}

// compute the residual
// n is the dimension
// M is the matrix
// the system is expressed as Ax = v; 
// input x here is the solution from the Jacobi solver
double *getResidual(double **M, double *v, double *x, int n){
    double *result = new double[n];
    double *vec = multiply_mv(M, n, n, x);
    for (int i = 0; i < n; i++){
        result[i] = v[i]-vec[i];
    }
    return result;
}

// V_cycle multigrid method for 2D problem
// A is the coefficient matrix
// v is the approximated solution
// f is the right hand side vector
// num is the number of nodes
// level is the current level of grid (from 1)
// pre and post are number of pre- and post-relaxation 
// num_level is the total number of grids 
double *V_cycle(double **A, double *v, double *f,\
 int row, int col, int level, int pre, int post, int num_level){
    // calculate the total number of nodes
    int num = row * col;

    // pre-relaxation
    A = five_stencil(row, col);
    v = Jac_p(num, pre, A, f, v);

    // compute the residual for the fine grid
    double *r = getResidual(A, f, v, num);

    // restriction 
    double *rc = restrict_2D(r, row, col);
    level = level + 1;
    int row_c = (row + 1) / 2;
    int col_c = (col + 1) / 2;
    int num_c = row_c * col_c;

    // inititalize the error vector
    double *e = new double[num_c];
    initialize_vec(e, num_c);

    // inititalize the coefficient matrix
    double **A_c = five_stencil(row_c, col_c);

    // stop recursion at coarsest grid, otherwise continue recursion
    if(level == num_level){
        //solve the error equation
        e = Jac_p(num_c, 1, A_c, rc, e);
        double *temp = getResidual(A_c, rc, e, num_c);
        double error = norm_max(temp, num_c);
        double conv = 1e-7;
        while (error > conv){
            e = Jac_p(num_c, 1, A_c, rc, e);
            temp = getResidual(A_c, rc, e, num_c);
            error = norm_max(temp, num_c);
        }
    } else{
        e = V_cycle(A_c, e, rc, row_c, col_c, level, pre, post, num_level);
    }

    // prolongation 
    double *ef = prolong_2D(e, row_c, col_c);

    // correct the approximated solution
    v = add_vv(v, ef, num);

    // post-relaxation
    v = Jac_p(num, post, A, f, v);

    return v;
}
