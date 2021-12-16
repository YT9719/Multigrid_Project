// header files for 2D problem 

#include <iostream>
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

// deallocate memory 
void deallocate_mat(double **A, int row, int col){
  for(int i = 0; i < row; i++){
    delete[] A[i];
  }
  delete[] A;
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
double **five_stencil(int row, int col, double dx){
  int i, j, d;
  d = row * col;
  double **m = allocate_mat(d, d);
  initialize_mat(m, d, d);
  for(int k = 0; k < d; k++){
      i = k/row; j = k-i*row;
      if(i>0)     {m[k][k-row] = 1.0/dx/dx;}
      if(i<col-1) {m[k][k+row] = 1.0/dx/dx;}
      if(j>0)     {m[k][k-1] = 1.0/dx/dx;}
      if(j<row-1) {m[k][k+1] = 1.0/dx/dx;}
      m[k][k] = -4.0/dx/dx;
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

// vector correction
// v and e are approximated solution and error
// n is the size of the vector
void correct(double *v, double *e, int n){
  for(int i = 0; i < n; i++){
    v[i] = v[i] + e[i];
  }
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
double *restrict_2D(double *rf, int row, int col,\
double *rc, int row_c, int col_c){
  int num = row * col; // number of nodes of fine grid 
  int num_c = row_c * col_c; // number of nodes of coarse grid
  double **R = allocate_mat(num_c, num);
  initialize_mat(R, num_c, num);
  int m, n;
  for(int i = 0; i < row_c; i++){
      for(int j = 0; j < col_c; j++){
            m = i * col_c + j;
            n = i * 2 * col + j * 2;
            R[m][n] = 1.0/16;
            R[m][n+1] = 2.0/16;
            R[m][n+2] = 1.0/16;
            R[m][n+col] = 2.0/16;
            R[m][n+col+1] = 4.0/16;
            R[m][n+col+2] = 2.0/16;
            R[m][n+col*2] = 1.0/16;
            R[m][n+col*2+1] = 2.0/16;
            R[m][n+col*2+2] = 1.0/16;
        }
    }
    rc = multiply_mv(R, num_c, num, rf);
    deallocate_mat(R, num_c, num);
    return rc;
}

// Prolongation module（2D）
// Full-weighting prolongation
// *ec is the error on the coarse grid
// row_c and col_c are number of rows and columns of coarse grid
// *ef is the error on the fine grid
double *prolong_2D(double *ec, int row_c, int col_c,\
double *ef, int row, int col){
    int num_c, num;
    num_c = row_c * col_c; // number of nodes of coarse grid
    num = row * col; // number of nodes of fine grid 
    double **P = allocate_mat(num, num_c);
    initialize_mat(P, num, num_c);
    int m, n;
    for(int i = 0; i < row_c; i++){
        for(int j = 0; j < col_c; j++){
            m = i * col_c + j;
            n = i * 2 * col + j * 2;
            P[n][m] = 1.0/4;
            P[n+1][m] = 2.0/4;
            P[n+2][m] = 1.0/4;
            P[n+col][m] = 2.0/4;
            P[n+col+1][m] = 4.0/4;
            P[n+col+2][m] = 2.0/4;
            P[n+col*2][m] = 1.0/4;
            P[n+col*2+1][m] = 2.0/4;
            P[n+col*2+2][m] = 1.0/4;
        }
    }
    ef = multiply_mv(P, num, num_c, ec);
    deallocate_mat(P, num, num_c);
    return ef;
}

// Gauss-Seidel solver
// n is the szie of vector
// **M is the coefficient matrix
// *v is the right hand side vector
// *x is the initial guess of solution 
double *GS(int n,int iter, double **M, double *v, double *x) {
    double sigma = 0;
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

// compute the residual
// n is the dimension
// A is the matrix
// r = f - Av
void residual(double *r, double **A, \
double *f, double *v, int n){
  double *temp = multiply_mv(A, n, n, v);
  for (int i = 0; i < n; i++){
    r[i] = f[i]-temp[i];
  }
  delete[] temp;
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
    
    // define necessary values and pointers
    int num = row * col;
    int row_c = (row-1)/2;
    int col_c = (col-1)/2;
    int num_c = row_c*col_c;
    double dx = 1.0/(row + 1);
    double dxc = dx*2;
    double *r = new double[num]; 
    double *rc = new double[num_c]; 
    double *temp = new double[num_c] ; 
    double *e = new double[num_c];
    double *ef = new double[num];

    // pre-relaxation
    A = five_stencil(row, col, dx);
    v = GS(num, pre, A, f, v);

    // compute the residual for the fine grid
    residual(r, A, f, v, num);

    // restriction 
    rc = restrict_2D(r, row, col, rc, row_c, col_c);
    level = level + 1;

    // inititalize the error vector
    initialize_vec(e, num_c);

    // inititalize the coefficient matrix
    double **A_c = five_stencil(row_c, col_c, dxc);

    // stop recursion at coarsest grid, otherwise continue recursion
    if(level == num_level){
        //solve the error equation
        e = GS(num_c, 1, A_c, rc, e);
        residual(temp, A_c, rc, e, num_c);
        double error = norm_max(temp, num_c);
        double conv = 1e-7;
        while (error > conv){
            e = GS(num_c, 1, A_c, rc, e);
            residual(temp, A_c, rc, e, num_c);
            error = norm_max(temp, num_c);
        }
    } else{
        e = V_cycle(A_c, e, rc, row_c, col_c, level, pre, post, num_level);
    }

    // prolongation 
    ef = prolong_2D(e, row_c, col_c, ef, row, col);

    // correct the approximated solution
    correct(v, ef, num);

    // post-relaxation
    v = GS(num, post, A, f, v);

    // free memory
    delete[] r, rc, e, ef, temp;

    return v;
}
