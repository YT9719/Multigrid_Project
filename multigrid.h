// include all modules for the multigrid method

#include <iostream>
#include <iomanip> // include the setprecision funtion
#include <cmath> // include math functions
#include <algorithm> // include max_element function
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

// generate coefficient matrix for three-point stencil
// n is the size of the square matrix 
double **three_stencil(int n){
  double **m = allocate_mat(n, n);
  initialize_mat(m, n, n);
  for(int i = 0; i < n; i++){
    if(i == 0){
      m[i][i] = 1; m[i][i+1] = -0.5;
    }
    else if(i == n - 1){
      m[i][i] = 1; m[i][i-1] = -0.5; 
    }
    else{
      m[i][i] = 1; m[i][i+1] = -0.5; m[i][i-1] = -0.5;
    }
  }
  return m;
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
    for(int i = 0; i < row; i++){
        double sum = 0;
        for(int j = 0; j < col; j++){
            sum = sum + m[i][j] * v[j];
        }
        b[i] = sum;
    }
    return b;
}

// vector subtraction
// *v1 & *v2 are two vectors
// n is the size of the vector
double *subtract_vv(double *v1, double *v2, int n){
  double *v3 = new double[n];
  for(int i = 0; i < n; i++){
    v3[i] = v1[i] - v2[i];
  }
  return v3;
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

// print vector
// *v is a vector
// n is the size of the vector
void print_v(double *v, int n){
    for(int i = 0; i < n; i++){
        cout<<fixed<<setprecision(5)<<v[i]<<endl;
    }
}

// print matrix
// **m is a vector
// row and col are number of rows and columns
void print_m(double **m, int row, int col){
    for(int i = 0; i < row; i++){
      for(int j = 0; j < col; j++ ){
        cout<<fixed<<setprecision(3)<<m[i][j]<<"  ";
      }
      cout<<endl;
    }
}

// Restriction module（1D）
// Full-weighting restriction
// *rf is the residual on the fine grid
// n is the szie of residual 
// *rc is the residual on the coarse grid
double *restrict(double *rf, int n){
    int l1 = n; // fine grid
    int l2 = (n + 1) / 2; // coarse grid
    double *rc = new double[l2];
    double **R = allocate_mat(l2, l1);
    initialize_mat(R, l2, l1);
    //for(int i = 0; i < l2; i++){R[i][i*2] = 1;} // same value
    for(int i = 0; i < l2; i++){
      if(i == 0){
        R[i][i*2] = 0.5;
        R[i][i*2+1] = 0.25;
      }
      else if(i == l2 - 1){
        R[i][2*i] = 0.5;
        R[i][2*i-1] = 0.25;
      }
      else{
        R[i][2*i-1] = 0.25;
        R[i][2*i] = 0.5;
        R[i][2*i+1] = 0.25;
      }
    }
    rc = multiply_mv(R, l2, l1, rf);
    delete[] R;
    return rc;
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


// Prolongation module
double *prolong(double *ec, int n){
    int l2 = n; // coarse grid
    int l1 = 2*n-1; // fine grid
    double *ef = new double [l1];
    double **P = allocate_mat(l1, l2);
    initialize_mat(P, l1, l2);
    // generate prolongation matrix
    for(int i = 0; i < l1; i++){
        if(i % 2 == 0){P[i][i/2] = 1;} // same value
        else{int j=(i - 1) / 2; 
        P[i][j] = 0.5;  // arithmetic average
        P[i][j+1] = 0.5;
        }
    }
    ef = multiply_mv(P, l1, l2, ec);
    delete[] P;
    return ef;
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



//Gauss Seidel Module
template <class T>
T* matMul(int n, T** M, T* v) {
  // n is the dimension of the linear system (matrix), M is the input matrix and v is the input vector
  int row, col;
  T* result = new T[n];
  for (row = 0; row < n; row++) {
    result[row] = 0;
    for (col = 0; col < n; col++) {
      result[row] += M[row][col]*v[col];
    }
  }
  return result;
}


// Gauss–Seidel solver
// n is the szie of vector
// **M is the coefficient matrix
// *v is the right hand side vector
// *x is the initial guess of solution 
template <class T>
T* GS(int n,int iter, T** M, T* v, T* x) {
  T sigma = 0;
  for (int k = 0; k < iter; k ++) {
    for (int i = 0; i < n; i ++) {
      sigma = 0;
      for (int j = 0; j < n; j ++) {
        if (j != i) {
          sigma += M[i][j]*x[j];
        }
      }
      x[i] = (v[i]-sigma)/M[i][i];
    }
  }
  return x;
}

//Get residual with vector x from GS procedure
template <class T>
T* getResidual(T** M,T* v, T* x, int n) {
 // n is the dimension
 // M is the matrix
 // the system is expressed as Ax = v; 
 // input x here is the solution from the GS procedure
 T* vec = new T[n];
 T* result = new T[n];
 vec = matMul<T>(n,M,x);
 for (int i = 0; i < n; i++){
   result[i] = v[i]-vec[i];
 }
 return result;
}
