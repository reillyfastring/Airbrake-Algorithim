#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

void choleskyDecomp(double* A, int n) {
    for (int i = 0; i < n; i++) {
        A[i*n + i] = sqrt(A[i*n + i] - dotProduct(i, &A[i*n], 1, &A[i*n], 1));
        for (int k = i + 1; k < n; k++) {
            A[k*n + i] = (A[k*n + i] - dotProduct(i, &A[k*n], 1, &A[i*n], 1))/ A[i*n + i];
        }
        for (int k = i + 1; k < n; k++) {
            A[i*n + k] = 0.0;
        }
    }
}

void scaleMatrix(double* A, int n, double scalar) {
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            A[i*n + k] = A[i*n + k]*scalar;
        }
    }
}

double dotProduct(const int N, const double *X, const int incX, const double *Y, const int incY) {
    double result = 0.0;
    int ix = 0;
    int iy = 0;

    for(int i = 0; i < N; i++) {
        result += X[ix] * Y[iy];
        ix += incX;
        iy += incY;
    }

    return result;
}
