/*
Basically, the prolongation(also called interpolation) 
transfers the data on a coarse grid to a fine grid. 
The simple way is to do a linear interplation to obtain
the data on the fine grid. 

The restriction refers to transfer the data on a find 
grid to a coarse grid. The simple way is to give the 
values on the fine grid to the coarse grid. 

Right now, the code focuses on a 1D problem.
*/
#define _USE_MATH_DEFINES //to include the mathematical constants
#include <iostream>
#include <iomanip> // include the setprecision funtion
#include <cmath> //include math functions

using namespace std;

// allocate a matrix 
double **allocate_mat(int row, int col){
    double** m = new double*[row];
    for(int i=0; i<row; i++){
        m[i] = new double[col];
    }
    return m;
}

// initialize a matrix
void initialize_mat(double **m, int row, int col)
{
    for(auto i = 0; i < row; ++i) {
        for(auto j = 0; j < col; ++j) {
            m[i][j] = 0;
        }
    }
}

// matrix and vector multiplication
double *multiply_mv(double **m, int row, int col, double *v){
    double *b = new double[row];
    for(int i=0; i<row; i++){
        double sum=0;
        for(int j=0; j<col; j++){
            sum = sum + m[i][j]*v[j];
        }
        b[i] = sum;
    }
    return b;
}

// print vector
void print_v(double* v, int n){
    for(int i=0; i<n; i++){
        cout<<fixed<<setprecision(5)<<v[i]<<endl;
    }
}

int main() {
    int n=4;
    double *x1 = new double[n]; //coarse grid
    double *x2 = new double[2*n-1]; // fine grid
    double *x3 = new double[n]; //back to the coarse grid

    // assign values for the coarse grid
    for(int i=0; i<n; i++){
        x1[i] = sin(i*M_PI/n);
    }

    // allocate prolongation matrix
    double** P = allocate_mat(2*n-1, n);
    initialize_mat(P, 2*n-1, n);

    // generate prolongation matrix
    for(int i=0; i<2*n-1; i++){
        if(i%2 == 0){P[i][i/2]=1;} // same value
        else{int j=(i-1)/2; P[i][j]=0.5; P[i][j+1]=0.5;} // arithmetic average
    }

    // print x1 vector
    cout<<"x1:"<<endl;
    print_v(x1,n);

    // apply prolongation matrix to the coarse grid
    x2 = multiply_mv(P, 2*n-1, n, x1);

    // print x2 vector
    cout<<"x2:"<<endl;
    print_v(x2,2*n-1);

    // allocate restriction matrix
    double** R = allocate_mat(n, 2*n-1);
    initialize_mat(R, n, 2*n-1);

    //generate restriction matrix
    for(int i=0; i<n; i++){
        R[i][i*2]=1; // same value
    }

    // apply restriction matrix to the coarse grid
    x3 = multiply_mv(R, n, 2*n-1, x2);

    // print x2 vector
    cout<<"x3:"<<endl;
    print_v(x3,n);

    return 1;
}

// Next step: apply prolongation and restriction to 2D problem

