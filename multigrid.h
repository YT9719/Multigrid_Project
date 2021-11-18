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
    for(auto i = 0; i < row; i++) {
        for(auto j = 0; j < col; j++) {
            m[i][j] = 0;
        }
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

// maximum norm of a vector
// *v is a vector
// n is the size of the vector
double norm_max(double *v, int n){
  for(int i = 0; i < n; i++){
    v[i] = abs(v[i]);
  }
  double max = *std::max_element(v, v+n);
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

// Restriction module 
// *rf is the residual on the fine grid
// n is the szie of residual 
// *rc is the residual on the coarse grid
double *restrict(double *rf, int n){ //n is the
    int l1 = n; // fine grid
    int l2 = (n + 1) / 2; // coarse grid
    double *rc = new double[l2];
    double **R = allocate_mat(l2, l1);
    initialize_mat(R, l2, l1);
    for(int i = 0; i < l2; i++){R[i][i*2] = 1;} // same value
    rc = multiply_mv(R, l2, l1, rf);
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
