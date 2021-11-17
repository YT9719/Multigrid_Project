// include all modules for the multigrid method

#include <iostream>
#include <iomanip> // include the setprecision funtion
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

// print vector
// *v is a vector
// n is the size of the vector
void print_v(double* v, int n){
    for(int i = 0; i < n; i++){
        cout<<fixed<<setprecision(5)<<v[i]<<endl;
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
T* GS(int n,int iter, T** M, T* v) {
  T* x = new T[n];
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
 T* vec = new T[n];
 T* result = new T[n];
 vec = matMul<T>(n,M,x);
 for (int i = 0; i < n; i++){
   result[i] = v[i]-vec[i];
 }
 return result
}
