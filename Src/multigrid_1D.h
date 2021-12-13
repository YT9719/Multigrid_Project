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

// generate coefficient matrix for three-point stencil
// n is the size of the square matrix 
double **three_stencil(int n, double dx){
  double **m = allocate_mat(n, n);
  initialize_mat(m, n, n);
  for(int i = 0; i < n; i++){
        if(i == 0){
            m[i][i] = -2.0/dx/dx; m[i][i+1] = 1.0/dx/dx;
        }
        else if(i == n - 1){
            m[i][i] = -2.0/dx/dx; m[i][i-1] = 1.0/dx/dx; 
        }
        else{
            m[i][i] = -2.0/dx/dx; m[i][i+1] = 1.0/dx/dx; m[i][i-1] = 1.0/dx/dx;
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

// Restriction module（1D）
// Full-weighting restriction
// *rf and *rc are the residual on the fine and coarse grid
// l1 and l2 are the number of nodes on the fine and coarse grid
double *restrict(double *rf, int l1, double *rc, int l2){
    // allocate and initialize restriction matrix
    double **R = allocate_mat(l2, l1);
    initialize_mat(R, l2, l1);
    // generate restriction matrix
    for(int i = 0; i < l2; i++){      
        R[i][2*i] = 0.25;
        R[i][2*i+1] = 0.5;
        R[i][2*i+2] = 0.25;
    }
    // apply the restriction matrix 
    rc = multiply_mv(R, l2, l1, rf);
    deallocate_mat(R, l2, l1);
    return rc;
}

// Prolongation module
// Full-weighting restriction
// *ef and *ec are the error on the fine and coarse grid
// l1 and l2 are the number of nodes on the fine and coarse grid
double *prolong(double *ec, int l2, double *ef, int l1){
    // allocate and initialize the prolongation matrix
    double **P = allocate_mat(l1, l2);
    initialize_mat(P, l1, l2);
    // generate prolongation matrix
    for(int i = 0; i < l2; i++){
        P[2*i][i]   = 0.5;
        P[2*i+1][i] = 1.0;
        P[2*i+2][i] = 0.5;
    }
    // apply prolongation matrix
    ef = multiply_mv(P, l1, l2, ec);
    deallocate_mat(P, l1, l2);
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

// V_cycle multigrid method for 1D problem
// A is the coefficient matrix
// v is the approximated solution
// f is the right hand side vector
// num is the number of nodes
// level is the current level of grid (from 1)
// pre and post are number of pre- and post-relaxation 
// num_level is the total number of grids 
double *V_cycle(double **A, double *v, double *f,\
int num, int level, int pre, int post, int num_level){
    
    // define necessary values and pointers
    int num_c = (num - 1) / 2; // number of nodes on the coarse grid 
    double dx = 1.0 / (num + 1); // grid spacing on the fine grid
    double dxc = dx * 2; // grid spacing on the coarse grid
    double *r = new double[num]; // residual pointer on the fine grid
    double *rc = new double[num_c]; // residual pointer on the coarse grid
    double *temp = new double[num_c] ; // temporary vector
    double *e = new double[num_c]; // error on the coarse grid
    double *ef = new double[num]; // error on the fine grid

    // pre-relaxation
    A = three_stencil(num, dx);
    v = GS(num, pre, A, f, v);

    // compute the residual for the fine grid
    residual(r, A, f, v, num);
    
    // restriction 
    rc = restrict(r, num, rc, num_c);
    level = level + 1;
    
    // inititalize the error vector
    initialize_vec(e, num_c);

    // inititalize the coefficient matrix
    double **A_c = three_stencil(num_c, dxc);

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
        e = V_cycle(A_c, e, rc, num_c, level, pre, post, num_level);
    }

    // prolongation 
    ef = prolong(e, num_c, ef, num);
    level = level - 1;

    // correct the approximated solution
    correct(v, ef, num);

    // post-relaxation
    v = GS(num, post, A, f, v);

    // free memory
    delete[] r, rc, e, ef, temp;

    return v;
}

// F-cycle multigrid code
double *F_cycle(double **A, double *v, double *f, \
int num, int level, int pre, int post, int num_level) {
    
    // define necessary values and pointers
    int num_c = (num - 1) / 2; // number of nodes on the coarse grid 
    double dx = 1.0 / (num + 1); // grid spacing on the fine grid
    double dxc = dx * 2; // grid spacing on the coarse grid
    double *r = new double[num]; // residual pointer on the fine grid
    double *rc = new double[num_c]; // residual pointer on the coarse grid
    double *temp = new double[num_c] ; // temporary vector
    double *e = new double[num_c]; // error on the coarse grid
    double *ef = new double[num]; // error on the fine grid

    // pre-smoothing
    A = three_stencil(num, dx);
    v = GS(num,pre,A,f,v);

    // compute the residual for the fine grid
    residual(r, A, f, v, num);

    // restriction
    rc = restrict(r, num, rc, num_c);
    level = level + 1;

    // initialize error vector
    initialize_vec(e,num_c);

    double **A_c = three_stencil(num_c, dxc);

    if (level == num_level) {
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
    }
    else {
        e = F_cycle(A_c,e,rc,num_c,level,pre,post,num_level);
    }

    ef = prolong(e, num_c, ef, num);
    level = level - 1;
    correct(v, ef, num);

    //re-smoothing
    v = GS(num, pre, A, f, v);
    residual(r, A, f, v, num);
  
    //restriction
    rc = restrict(r, num, rc, num_c);
    level = level + 1;

    // initialize error vector
    initialize_vec(e,num_c);

    if (level == num_level) {
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
    }
    else {
        e = F_cycle(A_c,e,rc,num_c,level,pre,post,num_level);
    }

    //prolongation
    ef = prolong(e, num_c, ef, num);
    level = level - 1;

    correct(v, ef, num);

    // free memory
    delete[] r, rc, e, ef, temp;

   // post relaxation
    v = GS(num,post,A,f,v);
    return v;
}