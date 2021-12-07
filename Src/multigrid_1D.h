// header files for 1D problem 

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

// Restriction module（1D）
// Full-weighting restriction
// *rf is the residual on the fine grid
// n is the szie of residual vector
// *rc is the residual on the coarse grid
double *restrict(double *rf, int n){

    // define the size of fine and coarse grid
    int l1 = n; // fine grid
    int l2 = (n + 1) / 2; // coarse grid

    // allocate and initialize residual vector and restriction matrix
    double *rc = new double[l2]; 
    double **R = allocate_mat(l2, l1);
    initialize_mat(R, l2, l1);

    // generate restriction matrix
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

    // apply the restriction matrix 
    rc = multiply_mv(R, l2, l1, rf);
    delete[] R;
    return rc;
}

// Prolongation module
// Full-weighting restriction
// *ec is the error on the coarse grid
// n is the szie of error vector 
// *ef is the the error on the fine grid
double *prolong(double *ec, int n){

    // define the size of coarse and fine grid
    int l2 = n; // coarse grid
    int l1 = 2*n-1; // fine grid

    // allocate and initialize error vector and prolongation matrix
    double *ef = new double [l1];
    double **P = allocate_mat(l1, l2);
    initialize_mat(P, l1, l2);

    // generate prolongation matrix
    for(int i = 0; i < l1; i++){
        if(i % 2 == 0){P[i][i/2] = 1;} 
        else{int j=(i - 1) / 2; 
        P[i][j] = 0.5;  
        P[i][j+1] = 0.5;
        }
    }

    // apply prolongation matrix
    ef = multiply_mv(P, l1, l2, ec);
    delete[] P;
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
// M is the matrix
// the system is expressed as Ax = v; 
// input x here is the solution from the GS solver
double *getResidual(double **M, double *v, double *x, int n){
    double *result = new double[n];
    double *vec = multiply_mv(M, n, n, x);
    for (int i = 0; i < n; i++){
        result[i] = v[i]-vec[i];
    }
    return result;
}

// A is the coefficient matrix
// v is the approximated solution
// f is the right hand side vector
// num is the number of nodes
// level is the current level of grid (from 1)
// pre and post are number of pre- and post-relaxation 
// num_level is the total number of grids 
double *V_cycle(double **A, double *v, double *f, \
int num, int level, int pre, int post, int num_level){
    
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
        e = V_cycle(A_c, e, rc, num_c, level, pre, post, num_level);
    }

    // prolongation 
    double *ef = prolong(e, num_c);

    // correct the approximated solution
    v = add_vv(v, ef, num);

    // post-relaxation
    v = GS(num, post, A, f, v);

    return v;
}


