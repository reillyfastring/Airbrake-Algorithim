#include "UKF_Operations.c"
#include "ProcessModel.c"
#define SIZE 9

int main() {
    //initialize state vector: 
    double StateVector[SIZE] = {0};
    double StateMean[SIZE] = {0};
    
    //initialize sensor vector:
    //gryo, accel, and mag can poll at 100Hz, Gps can poll at 1.15Hz
    double sensorVector[SIZE] = {0};

    //initializing covariance matrix
    double Covariance[SIZE*SIZE] = {0};
    for(int i = 0; i < SIZE; i++) {
        Covariance[i*SIZE + i] = 1; 
    }
    double PredictedCovariance[SIZE*SIZE] = {0};
    for(int i = 0; i < SIZE; i++) {
        Covariance[i*SIZE + i] = 1; 
    }


    //sigma points
    double sigma[(2*SIZE+1)*SIZE];
    double sigmaPredicted[(2*SIZE+1)*SIZE];
    double weights[(2*SIZE+1)*2];

    //Sigma point tuning parameters
    double a = 1*pow(10,-3); //spread of sigma points around mean
    double b = 2; //uses prior knowledge about distrubution of state Gaussian is usually 2
    double k = -9; //scaling parameter
    
    //Sigma Point calculation
    double y = pow(a, 2)*(SIZE+k)-SIZE;

    choleskyDecomp(Covariance, SIZE);
    scaleMatrix(Covariance, SIZE, sqrt(SIZE+y));

    //sets all sigma points to the state vector
    for (int i = 0; i < (2*SIZE+1); i++) {
        for (int k = 0; k < SIZE; k++) {
            sigma[i*SIZE + k] = StateVector[k];
        }
    }
    //adds covariance columns for half of sigma points
    for (int i = 1; i < (SIZE+1); i++) {
        for (int k = 0; k < SIZE; k++) {
            sigma[i*SIZE + k] = sigma[i*SIZE + k] + Covariance[k*SIZE+(i-1)];
        }
    }
    //subtracts covariance columns for other half of sigma points
    for (int i = (SIZE+1); i < (2*SIZE+1); i++) {
        for (int k = 0; k < SIZE; k++) {
            sigma[i*SIZE + k] = sigma[i*SIZE + k] - Covariance[k*SIZE +(i-13)];
        }
    }
    //calculates weights
    weights[0] = y/(SIZE+y);
    weights[1] = weights[0] + (1-pow(a, 2)+b);
    for (int i = 2; i < (2*SIZE+1); i++) {
        weights[i] = 1/(2*(SIZE+y));
    }

    //Propogate Sigma Points


}
