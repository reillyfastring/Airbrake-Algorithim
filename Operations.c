#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

double randn() {
    double u1 = (rand() + 1.0) / (RAND_MAX + 1.0); 
    double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

void drag(double theta, double velocity) {
//function that calculates drag using model points, surface area, and current mach speed.
}

double rho(double height) {
    return (pow(2116*((59+459.7-.00356*(height*3.281))/518.6), 5.256)/(1718*(59+459.7-.00356*(height*3.281)))*515.4);
}

double surfaceA(angleOfDeployment) {
    surfaceA = 0.0182414692475+angleofDeployment*0.00350967*9;
}
