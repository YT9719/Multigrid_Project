/*
Basically, the prolongation(also called interpolation) 
transfers the data on a coarse grid to a fine grid. 
The simple way is to do a linear interplation to obtain
the data on a fine grid. 
*/
#define _USE_MATH_DEFINES //to include the mathematical constants
#include <iostream>
#include <cmath> //include math functions

using namespace std;

// allocate a matrix 
double **allocate_mat(int row, int col){
    double** M = new double*[row];
    for(int i=0; i<row; i++){
        M[i] = new double[col];
    }
    return M;
}

// initialize matrix
void initialize_mat(double **m, int row, int col)
{
    for(auto i = 0; i < row; ++i) {
        for(auto j = 0; j < col; ++j) {
            m[i][j] = 0;
        }
    }
}

int main() {
    int n=4;
    double *x1 = new double[n]; //coarse grid
    double *x2 = new double[2*n-1]; // fine grid

    // assign values for the coarse grid
    for(int i=0; i<n; i++){
        x1[i] = sin(i*M_PI/n);
    }

    // allocate interpolation matrix
    double** I_ctof = allocate_mat(2*n-1, n);
    initialize_mat(I_ctof, 2*n-1, n);

    // generate interpolation matrix
    for(int i=0; i<2*n-1; i++){
        if(i%2 == 0){I_ctof[i][i/2]=1;} // same value
        else{int j=(i-1)/2; I_ctof[i][j]=0.5; I_ctof[i][j+1]=0.5;} // arithmetic average
    }

    // print interpolation matrix
    for(int i=0; i<2*n-1; i++){
        for(int j=0; j<n; j++){
            cout<<I_ctof[i][j]<<"   ";
        }
        cout<<endl;
    }

    return 1;
}

