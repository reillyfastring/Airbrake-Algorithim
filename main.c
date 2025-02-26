#include "ukf.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// You might define or retrieve these from somewhere else:
static const double t = 0.01;             // Timestep
static double DeployAngle = 0.0;          // Some initial angle if needed

int main(void)
{
    // ============= Initialization =============
    
    // State vector: 
    // [ dragZ, velZ, posZ, thetaX, thetaY, thetaZ,
    //   accelbias, gyrobias, gpsbias,
    //   accelX, accelY, accelZ,
    //   omegaX, omegaY, omegaZ,
    //   gpsX, gpsY, gpsZ ]
    double StateVector[SIZE] = {0.0};

    double StateMean[SIZE] = {0.0};
    double MeasurementMean[SIZE] = {0.0};   // If your measurement dimension = SIZE/2,
                                              // you can store it in the first half 
    double previousAngle = 0.0;
    double Angle = 0.0;
    
    // Sensor vector: e.g. [accelX, accelY, accelZ, gyroX, gyroY, gyroZ, gpsX, gpsY, gpsZ]
    // or likewise if it’s 9 elements. Here we keep SIZE=18 for consistency.
    double sensorVector[SIZE] = {0.0};

    // Covariance matrices
    double Covariance[SIZE*SIZE] = {0.0};
    double PredictedCovariance[SIZE*SIZE] = {0.0};
    double PredictedMeasurementCovariance[SIZE*SIZE] = {0.0};
    double crossCovariance[(SIZE)*(SIZE/2)] = {0.0};  // For state x measurement

    // Process & Measurement noise
    double ProcessNoise[SIZE*SIZE] = {0.0};
    double MeasurementNoise[(SIZE/2)*(SIZE/2)] = {0.0};

    // Kalman Gain: (SIZE x SIZE/2)
    double kalmanGain[(SIZE)*(SIZE/2)] = {0.0};

    // Sigma points
    double sigma[(2*SIZE + 1) * SIZE] = {0.0};
    double measurementMatrix[(2*SIZE + 1)*(SIZE/2)] = {0.0};

    // Weight vector
    // Some references show that you need two arrays: Wm (mean) & Wc (cov), but 
    // here you store them in one array or in two. We keep it simple:
    double weights[2*SIZE + 1] = {0.0};
    
    // ============= Setup initial covariance & noise =============

    // Initialize main covariance to identity
    for(int i = 0; i < SIZE; i++){
        Covariance[i*SIZE + i] = 1.0;
        PredictedCovariance[i*SIZE + i] = 1.0;
        PredictedMeasurementCovariance[i*SIZE + i] = 1.0;
    }

    // Process noise = identity
    for(int i = 0; i < SIZE; i++){
        ProcessNoise[i*SIZE + i] = 1.0;
    }

    // Measurement noise (9x9 if SIZE/2=9)
    for(int i = 0; i < SIZE/2; i++){
        MeasurementNoise[i*(SIZE/2) + i] = 1.0;
    }

    // ============= Sigma point tuning parameters =============
    double a = 1e-3;  // spread of sigma points
    double b = 2.0;   // knowledge of distribution (2 for Gaussian)
    double k = -9.0;  // scaling
    double y = pow(a, 2) * (SIZE + k) - SIZE;  // \lambda = alpha^2 (n + k) - n

    // ============= Cholesky-based scaling of covariance =============
    choleskyDecomp(Covariance, SIZE);
    scaleMatrix(Covariance, SIZE, sqrt(SIZE + y));

    // ============= Generate sigma points from StateVector =============
    // Each row i => sigma point i
    // sigma[i*SIZE + col], col in [0..SIZE-1]
    for(int i = 0; i < (2*SIZE + 1); i++) {
        for(int c = 0; c < SIZE; c++) {
            sigma[i*SIZE + c] = StateVector[c];
        }
    }
    // Add columns for first half
    for(int i = 1; i <= SIZE; i++) {
        for(int c = 0; c < SIZE; c++) {
            sigma[i*SIZE + c] += Covariance[c*SIZE + (i-1)];
        }
    }
    // Subtract columns for second half
    for(int i = SIZE+1; i < (2*SIZE + 1); i++) {
        for(int c = 0; c < SIZE; c++) {
            // index (i-13) might be an offset for your scenario,
            // but be sure that's correct for your dimension. 
            sigma[i*SIZE + c] -= Covariance[c*SIZE + (i - (SIZE+1))];
        }
    }

    // ============= Weights =============
    weights[0] = y / (SIZE + y);
    // some references use separate mean/cov weights, but we'll keep it simple
    // The second entry is typically W_0^c = W_0^m + (1 - alpha^2 + beta)
    // or something similar. For a minimal fix, do:
    weights[1] = weights[0] + (1 - pow(a, 2) + b);

    for(int i = 2; i < (2*SIZE+1); i++){
        weights[i] = 1.0 / (2.0 * (SIZE + y));
    }

    // ============= Start Filter Steps =============

    // 1) ProcessModel
    ProcessModel(sigma, t, &DeployAngle, sensorVector);

    // 2) Compute Mean
    computeMean(sigma, weights, StateMean);

    // 3) Predict Covariance
    PredictCovariance(sigma, StateMean, PredictedCovariance, weights, ProcessNoise);

    // 4) Measurement
    MeasurementFunction(measurementMatrix, sigma);

    // 5) Compute Measurement Mean
    computeMean(measurementMatrix, weights, MeasurementMean);

    // 6) Predict Measurement Covariance
    PredictCovariance(measurementMatrix, MeasurementMean, PredictedMeasurementCovariance, weights, MeasurementNoise);
    choleskyDecomp(PredictedMeasurementCovariance, SIZE);

    // 7) Cross Covariance
    CrossCovariance(sigma, measurementMatrix, StateMean, MeasurementMean, crossCovariance, weights);

    // 8) Kalman Gain => K = crossCov * inv(PmeasCov)
    invertMatrixCholesky(PredictedMeasurementCovariance, PredictedMeasurementCovariance, SIZE);
    matrixMultiply(crossCovariance, PredictedMeasurementCovariance, kalmanGain, SIZE, SIZE/2, SIZE);
    matrixTranspose(PredictedMeasurementCovariance, PredictedMeasurementCovariance, SIZE);
    matrixMultiply(kalmanGain, PredictedMeasurementCovariance, kalmanGain, SIZE, SIZE/2, SIZE);

    // 9) Update state vector
    updateState(StateMean, kalmanGain, MeasurementMean, sensorVector, StateVector);

    // 10) update of covariance matrix
    cholUpdateMulti(covariance, measurementMatrix, SIZE, SIZE/2, -1);

    // 11) Suppose we compute next angle:
    previousAngle = Angle;
    Angle = PredictDeploymentAngle(StateVector, previousAngle);

    printf("Updated angle = %f\n", Angle);
    printf("Updated state example, posZ = %f\n", StateVector[2]);

    return 0;
}
