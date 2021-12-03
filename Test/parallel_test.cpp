#define _USE_MATH_DEFINES //to include the mathematical constants
#include <iostream>
#include <iomanip> // include the setprecision functions
#include <cmath> // include math functions
#include "../multigrid.h" // include all modules and functions
#include <omp.h> // include OpenMP

using namespace std;

#define MAX_THREADS 4

double *multiply(double **m, int row, int col, double *v){
    double *b = new double [row];

    #pragma omp parallel 
    {   
        #pragma omp for
        for(int i = 0; i < row; i++){
            double sum = 0;
            for(int j = 0; j < col; j++){
                sum = sum + m[i][j] * v[j];
            }
            b[i] = sum;
        }
    }

    return b;
}

int main(){

    int n = 10000;
    int num = n + 1;
    double **A = three_stencil(num);
    double *x  = new double[num];
    initialize_vec(x, num);
    double *f  = new double[num];
    
    double start_time = omp_get_wtime();

    f =  multiply(A, num, num, x);

    double run_time = omp_get_wtime() - start_time;

    cout<<run_time<<endl;
    //cout<<"f:"<<endl;
    //print_v(f, num);

    return 1;
}