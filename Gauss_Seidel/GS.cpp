#include <iostream>
using namespace std;

template <class T>
T* matMul(int n, T** M, T* v) {
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
/*

template <class T>
T* GS(int n, T** M, T* v) {

}
*/
int main() {
  int n = 2;
  double** M = new double*[n];
  
  for (int i = 0; i < n; i ++) {
    M[i] = new double[n];
  } 
  
  double *v = new double[n];

  for(int i = 0; i < n; i ++) {
    v[i] = i+n;
    for (int j = 0; j < n; j++) {
      M[i][j] = i+j+n;
    }
  }
  M[0][0] = 16;
  M[0][1] = 3;
  M[1][0] = 7;
  M[1][1] = -11;

  v[0] = 11;
  v[1] = 13;
  double* x = GS(n, 100, M, v); 
 
  double* result = matMul<double>(n,M,v);  

  for(int i = 0; i < n; i ++) {
    for (int j = 0; j < n; j++) {
      cout << M[i][j] << " ";
    }
    cout << endl;
  }
/*
  for (int i = 0; i < n; i ++) {
    cout << v[i] << "  "; 
    cout << result[i] << endl;
  }
*/
  for (int i = 0; i < n; i ++) {
    cout << x[i] << endl;
  }  
  for (int i = 0; i < n; i ++) {
    delete[] M[i];
  }
  delete[] M;
  delete[] v;
  delete[] result;
  delete[] x;
  return 0;
}
