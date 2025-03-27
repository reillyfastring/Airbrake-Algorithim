/* ukf.c
 * Implementation of Unscented Kalman Filter helper functions.
 */

 #include "ukf.h"
 #include <math.h>
 #include <stdlib.h>
 #include <stdio.h>
 
 // Optionally define PI if not defined
 #ifndef M_PI
 #define M_PI 3.14159265358979323846
 #endif
 
 //--------------------------------------------------//
 // Process Model
 //--------------------------------------------------//
 void processModel(double* sigma, double t, double* DeployAngle, double* sensorVector)
 {
     // dragZ(0), velocityZ(1), positionZ(2), thetaX(3), thetaY(4), thetaZ(5),
     // accelBias(6), gyroBias(7), gpsBias(8), sensor placeholders(9..17)
     for (int i = 0; i < (2*SIZE + 1); i++) {
         // Example: compute magnitude of accel from sensor vector
         double magA = sqrt(
             pow(sensorVector[1], 2) +
             pow(sensorVector[2], 2) +
             pow(sensorVector[3], 2)
         );
 
         // Prediction for dragZ
         sigma[i*SIZE + 0] += drag(*DeployAngle, sigma[i*SIZE + 1]);
 
         // velocityZ
         sigma[i*SIZE + 1] += sin(magA - sigma[i*SIZE + 0]) * t;
 
         // positionZ
         sigma[i*SIZE + 2] += sigma[i*SIZE + 1] * t;
 
         // thetaX, thetaY, thetaZ (using sensorVector as raw rates?)
         sigma[i*SIZE + 3] += sensorVector[3] * t;
         sigma[i*SIZE + 4] += sensorVector[4] * t;
         sigma[i*SIZE + 5] += sensorVector[5] * t;
 
         // accelerate / gyro / GPS bias random walk
         sigma[i*SIZE + 6] += AccC * randn();
         sigma[i*SIZE + 7] += GyroC * randn();
         sigma[i*SIZE + 8] += GPSC * randn();
 
         // store sensor readings in last 9 slots
         for (int s = 0; s < 9; s++) {
             sigma[i*SIZE + (9 + s)] = sensorVector[s];
         }
     }
 }
 
 //--------------------------------------------------//
 // Compute Mean of Sigma Points
 //--------------------------------------------------//
 void computeMean(double* sigma, double* Weights, double* stateMean)
 {
     // Zero out
     for (int j = 0; j < SIZE; j++) {
         stateMean[j] = 0.0;
     }
 
     // Weighted sum
     for (int i = 0; i < (2*SIZE + 1); i++) {
         double w = Weights[i];
         for (int j = 0; j < SIZE; j++) {
             stateMean[j] += w * sigma[i*SIZE + j];
         }
     }
 }
 
 //--------------------------------------------------//
 // Predict State Covariance
 //--------------------------------------------------//
 void PredictCovariance(double* sigma, double* stateMean,
                        double* PredictedCovariance, double* weights,
                        double* ProcessNoise)
 {
     // Start with process noise
     for (int idx = 0; idx < SIZE*SIZE; idx++) {
         PredictedCovariance[idx] = ProcessNoise[idx];
     }
 
     // Sum contributions
     for (int i = 0; i < (2*SIZE + 1); i++) {
         double w = weights[i];
         double difference[SIZE];
         for (int k = 0; k < SIZE; k++) {
             difference[k] = sigma[i*SIZE + k] - stateMean[k];
         }
         for (int r = 0; r < SIZE; r++) {
             for (int c = 0; c < SIZE; c++) {
                 PredictedCovariance[r*SIZE + c] += w * difference[r] * difference[c];
             }
         }
     }
 }
 
 //--------------------------------------------------//
 // Measurement Function
 //--------------------------------------------------//
 void MeasurementFunction(double* Measurement, double* sigma)
 {
     // Suppose measurement dimension = SIZE/2 = 9.
     // The last 9 indices of the state: sensor(9..17)
     // plus biases(6..8)? Here we add them as an example
     for (int i = 0; i < (2*SIZE + 1); i++) {
         // accelX = sigma[i*SIZE + 9] + accelBias(6)
         Measurement[i*(SIZE/2) + 0] = sigma[i*SIZE +  9] + sigma[i*SIZE + 6];
         Measurement[i*(SIZE/2) + 1] = sigma[i*SIZE + 10] + sigma[i*SIZE + 6];
         Measurement[i*(SIZE/2) + 2] = sigma[i*SIZE + 11] + sigma[i*SIZE + 6];
 
         // gyroX = sigma[i*SIZE + 12] + gyroBias(7)
         Measurement[i*(SIZE/2) + 3] = sigma[i*SIZE + 12] + sigma[i*SIZE + 7];
         Measurement[i*(SIZE/2) + 4] = sigma[i*SIZE + 13] + sigma[i*SIZE + 7];
         Measurement[i*(SIZE/2) + 5] = sigma[i*SIZE + 14] + sigma[i*SIZE + 7];
 
         // gpsX = sigma[i*SIZE + 15] + gpsBias(8)
         Measurement[i*(SIZE/2) + 6] = sigma[i*SIZE + 15] + sigma[i*SIZE + 8];
         Measurement[i*(SIZE/2) + 7] = sigma[i*SIZE + 16] + sigma[i*SIZE + 8];
         Measurement[i*(SIZE/2) + 8] = sigma[i*SIZE + 17] + sigma[i*SIZE + 8];
     }
 }
 
 //--------------------------------------------------//
 // Cross Covariance
 //--------------------------------------------------//
 void CrossCovariance(double* sigma, double* measurementMatrix,
                      double* stateMean, double* measurementMean,
                      double* crossCovariance, double* weights)
 {
     // crossCov is (SIZE x SIZE/2)
     for (int idx = 0; idx < SIZE*(SIZE/2); idx++) {
         crossCovariance[idx] = 0.0;
     }
 
     for (int i = 0; i < (2*SIZE + 1); i++) {
         double w = weights[i];
         for (int c = 0; c < SIZE; c++) {
             double diffState = sigma[i*SIZE + c] - stateMean[c];
             for (int d = 0; d < SIZE/2; d++) {
                 double diffMeas = measurementMatrix[i*(SIZE/2) + d] - measurementMean[d];
                 crossCovariance[c*(SIZE/2) + d] += w * diffState * diffMeas;
             }
         }
     }
 }
 
 //--------------------------------------------------//
 // Cholesky Update Multi
 //--------------------------------------------------//
 void cholUpdateMulti(double* Sx, double* U, int n, int m, int sign)
 {
     double temp[SIZE];
     for (int j = 0; j < m; j++) {
         for (int i = 0; i < n; i++) {
             temp[i] = U[i*m + j];
         }
         cholUpdate(Sx, temp, n, sign);
     }
 }
 
 //--------------------------------------------------//
 // Cholesky Rank-1 Update
 //--------------------------------------------------//
 void cholUpdate(double* S, double* x, int n, int sign)
 {
     // This modifies S (lower-tri) to reflect +/- x*x^T in S*S^T.
     for (int i = 0; i < n; i++) {
         double Sii = S[i*n + i];
         double xi  = x[i];
 
         double alpha = 0.0;
         if (sign > 0) {
             alpha = Sii*Sii + xi*xi;  // update
         } else {
             alpha = Sii*Sii - xi*xi;  // downdate
         }
 
         // We might clamp alpha if <= 0
         if (alpha < 1e-14) {
             alpha = 1e-14;
         }
 
         double r = sqrt(alpha);
         double c = r / Sii;
         double s = xi / Sii;
 
         S[i*n + i] = r;
 
         // Off-diagonal
         for (int j = i+1; j < n; j++) {
             double Sj_i = S[j*n + i];
             double beta = (sign > 0)
                           ? (Sj_i + s*x[j])
                           : (Sj_i - s*x[j]);
             S[j*n + i] = beta / c;
 
             // Next iteration modifies x[j]
             x[j] = c*x[j] - s*Sj_i;
         }
     }
 }
 
 //--------------------------------------------------//
 // Update State
 //--------------------------------------------------//
 void updateState(double* X0, double* K, double* Zpred,
                  double* Zinit, double* stateVector)
 {
     // measurement dimension = SIZE/2
     double error[SIZE/2];
     for (int i = 0; i < SIZE/2; i++) {
         error[i] = Zinit[i] - Zpred[i];
     }
 
     // Correction = K * error
     double correction[SIZE];
     matrixMultiply(K, error, correction, SIZE, SIZE/2, 1);
 
     // stateVector = X0 + correction
     for (int i = 0; i < SIZE; i++) {
         stateVector[i] = X0[i] + correction[i];
     }
 }
 
 //--------------------------------------------------//
 // Basic Linear Algebra Helpers
 //--------------------------------------------------//
 void choleskyDecomp(double* A, int n)
 {
     // in-place Cholesky: A -> lower tri
     for (int i = 0; i < n; i++) {
         double diag = A[i*n + i]
                       - dotProduct(i, &A[i*n], 1, &A[i*n], 1);
         if (diag < 1e-14) {
             diag = 1e-14;
         }
         A[i*n + i] = sqrt(diag);
 
         for (int k = i+1; k < n; k++) {
             double val = A[k*n + i]
                          - dotProduct(i, &A[k*n], 1, &A[i*n], 1);
             A[k*n + i] = val / A[i*n + i];
         }
         // zero out upper
         for (int k = i+1; k < n; k++) {
             A[i*n + k] = 0.0;
         }
     }
 }
 
 void scaleMatrix(double* A, int n, double scalar)
 {
     for (int i = 0; i < n; i++){
         for (int j = 0; j < n; j++){
             A[i*n + j] *= scalar;
         }
     }
 }
 
 void invertLowerTriangular(double* L, double* L_inv, int n)
 {
     // Invert lower-triangular L -> L_inv
     for (int i = 0; i < n; i++) {
         L_inv[i*n + i] = 1.0 / L[i*n + i];
         for (int j = 0; j < i; j++){
             double sum = 0.0;
             for (int k = j; k < i; k++){
                 sum += L[i*n + k] * L_inv[k*n + j];
             }
             L_inv[i*n + j] = -sum / L[i*n + i];
         }
         for (int j = i+1; j < n; j++){
             L_inv[i*n + j] = 0.0;
         }
     }
 }
 
 void multiplyTranspose(double* L_inv, double* A_inv, int n)
 {
     // A_inv = L_inv^T * L_inv
     for (int i = 0; i < n; i++){
         for (int j = 0; j < n; j++){
             double sum = 0.0;
             for (int k = 0; k < n; k++){
                 sum += L_inv[k*n + i] * L_inv[k*n + j];
             }
             A_inv[i*n + j] = sum;
         }
     }
 }
 
 void invertMatrixCholesky(double* A, double* A_inv, int n)
 {
     double* L     = (double*)malloc(n*n*sizeof(double));
     double* L_inv = (double*)malloc(n*n*sizeof(double));
 
     // copy A -> L
     for(int i = 0; i < n*n; i++){
         L[i] = A[i];
     }
 
     // Decompose
     choleskyDecomp(L, n);
 
     // Invert L
     invertLowerTriangular(L, L_inv, n);
 
     // A_inv = L_inv^T * L_inv
     multiplyTranspose(L_inv, A_inv, n);
 
     free(L);
     free(L_inv);
 }
 
 //--------------------------------------------------//
 // Dot Product
 //--------------------------------------------------//
 double dotProduct(const int N, const double* X, const int incX,
                   const double* Y, const int incY)
 {
     double result = 0.0;
     int ix = 0;
     int iy = 0;
     for (int i = 0; i < N; i++){
         result += X[ix] * Y[iy];
         ix += incX;
         iy += incY;
     }
     return result;
 }
 
 //--------------------------------------------------//
 // Matrix Multiply
 //--------------------------------------------------//
 void matrixMultiply(double* A, double* B, double* C,
                     int rowsA, int colsA, int colsB)
 {
     // C = A * B
     for (int i = 0; i < rowsA; i++){
         for (int j = 0; j < colsB; j++){
             double sum = 0.0;
             for (int k = 0; k < colsA; k++){
                 sum += A[i*colsA + k] * B[k*colsB + j];
             }
             C[i*colsB + j] = sum;
         }
     }
 }
 
 //--------------------------------------------------//
 // Matrix Transpose (square)
 //--------------------------------------------------//
 void matrixTranspose(double* A, double* A_T, int a)
 {
     for (int i = 0; i < a; i++){
         for (int j = 0; j < a; j++){
             A_T[j*a + i] = A[i*a + j];
         }
     }
 }
 
 //--------------------------------------------------//
 // randn() Implementation
 //--------------------------------------------------//
 double randn(void)
 {
     double u1 = (rand() + 1.0) / ((double)RAND_MAX + 1.0);
     double u2 = (rand() + 1.0) / ((double)RAND_MAX + 1.0);
     return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
 }
 
 //--------------------------------------------------//
 // Example drag
 //--------------------------------------------------//
 double drag(double theta, double velocity)
 {
     // dummy function returning a constant or some small model
    return(-0.06828*velocity/343*0.01+.243866*theta+.333907);
 }
 
 //--------------------------------------------------//
 // Example density model
 //--------------------------------------------------//
 double rho(double height)
 {
     // Possibly a mock
     double T = 15.04 - 0.00649 * height;
     double p = 101.29 * pow((T + 273.1) / 288.08, 5.256);
     double rho = p / (0.2869 * (T + 273.1));
     return rho;
 }
 
 //--------------------------------------------------//
 // Example surface area
 //--------------------------------------------------//
 double surfaceA(double angleOfDeployment)
 {
     // return 0.01824 + angleOfDeployment*someFactor;
     return 0.0182414692475 + angleOfDeployment*0.00350967*9.0;
 }
 
 //--------------------------------------------------//
 // PredictDeploymentAngle
 //--------------------------------------------------//
 double PredictDeploymentAngle(double* stateVector)
 {
     double low = 0.0;
     double high = 1.0;
     double mid = 0.0;
     double targetApogee = 0;
     if (stateVector[2] < 100) {
        targetApogee = 9144;
     }
     else {
        targetApogee = 9144 + (9144 - stateVector[3])/20;
     }
     // Keep track of how many times we call PredictApogee
     int predictCalls = 0;
     double predictedApogee = 0;
 
     // Continue until our angle interval is small enough
     while ((high - low) > 0.0001)
     {
         mid = (low + high) / 2.0;
         
         // Call PredictApogee to check how close we are to the target
         predictedApogee = PredictApogee(stateVector, mid);
         predictCalls++;
 
         // If we overshoot 9144, we shrink the upper bound; 
         // otherwise, we shrink the lower bound.
         if (predictedApogee > targetApogee)
         {
             low = mid;
         }
         else
         {
             high = mid;
         }
     }
 
     // Print how many times we called PredictApogee
     printf("PredictApogee was called %d times.\nThe final Predicted Apogee was %f\n", predictCalls, predictedApogee);
 
     // Return the midpoint of our final [low, high] interval as the "best" angle
     return (low + high) / 2.0;
 }
 
 
 double PredictApogee(double* stateVector, double deploymentAngle)
 {
     // Retrieve initial state variables
     double velocityZ = stateVector[2];
     double positionZ = stateVector[3];
     double thetaZ    = stateVector[6];
     
     // Initial timestep guess
     double dt = 0.3;
     int maxIters = 10000000;
     int iter = 0;
 
     while(velocityZ > 0.0 && iter < maxIters) {
         // k1
         if (stateVector[6]*180/M_PI < 20) {
            deploymentAngle = 0;
         }
         double k1_rho   = rho(positionZ);
         double k1_x     = velocityZ;
         double k1_v     = -grav
                           - (0.5 / MASS) * k1_rho * drag(deploymentAngle, velocityZ)
                           * surfaceA(deploymentAngle)
                           * pow(velocityZ, 2)
                           * (1.0 / sin(thetaZ));
         double k1_theta = grav * sin(thetaZ) * cos(thetaZ) / velocityZ;
 
         // k2
         double vk1    = velocityZ + 0.5 * dt * k1_v;
         double posk1  = positionZ + 0.5 * dt * k1_x;
         double thetaK1 = thetaZ + 0.5 * dt * k1_theta;
         
         double k2_rho   = rho(posk1);
         double k2_x     = vk1;
         double k2_v     = -grav
                           - (0.5 / MASS) * k2_rho * drag(deploymentAngle, vk1)
                           * surfaceA(deploymentAngle)
                           * pow(vk1, 2)
                           * (1.0 / sin(thetaK1));
         double k2_theta = grav * sin(thetaK1) * cos(thetaK1) / vk1;
 
         // k3
         double vk2    = velocityZ + 0.5 * dt * k2_v;
         double posk2  = positionZ + 0.5 * dt * k2_x;
         double thetaK2 = thetaZ + 0.5 * dt * k2_theta;
         
         double k3_rho   = rho(posk2);
         double k3_x     = vk2;
         double k3_v     = -grav
                           - (0.5 / MASS) * k3_rho * drag(deploymentAngle, vk2)
                           * surfaceA(deploymentAngle)
                           * pow(vk2, 2)
                           * (1.0 / sin(thetaK2));
         double k3_theta = grav * sin(thetaK2) * cos(thetaK2) / vk2;
 
         // k4
         double vk3    = velocityZ + dt * k3_v;
         double posk3  = positionZ + dt * k3_x;
         double thetaK3 = thetaZ + dt * k3_theta;
         
         double k4_rho   = rho(posk3);
         double k4_x     = vk3;
         double k4_v     = -grav
                           - (0.5 / MASS) * k4_rho * drag(deploymentAngle, vk3)
                           * surfaceA(deploymentAngle)
                           * pow(vk3, 2)
                           * (1.0 / sin(thetaK3));
         double k4_theta = grav * sin(thetaK3) * cos(thetaK3) / vk3;
 
         // Update the state using 4th-order Runge-Kutta
         positionZ += (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x) * dt / 6.0;
         velocityZ += (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v) * dt / 6.0;
         thetaZ    -= (k1_theta + 2.0*k2_theta + 2.0*k3_theta + k4_theta) * dt / 6.0;
 
         // Clamp dt_candidate between 0.01 and 1.5 seconds.
         if (velocityZ > 343) {
             dt = .8;
         }
         if (velocityZ < 50) {
             dt = 0.2;
         }
         
         iter++;
         // Optionally print out iteration details for debugging:
         // printf("Iteration %d: posZ = %f, velocityZ = %f, dt = %f\n", iter, positionZ, velocityZ, dt);
     }
     
     //printf("Iterations: %d\n", iter);
     return positionZ;
 }

 double* PredictApogeeSS(double* stateVector, double deploymentAngle, double dt)
 {
     // Retrieve initial state variables
     double velocityZ = stateVector[2];
     double positionZ = stateVector[3];
     double thetaZ    = stateVector[6];
     
         // k1
         double k1_rho   = rho(positionZ);
         double k1_x     = velocityZ;
         double k1_v     = -grav
                           - (0.5 / MASS) * k1_rho * drag(deploymentAngle, velocityZ)
                           * surfaceA(deploymentAngle)
                           * pow(velocityZ, 2)
                           * (1.0 / sin(thetaZ));
         double k1_theta = grav * sin(thetaZ) * cos(thetaZ) / velocityZ;
 
         // k2
         double vk1    = velocityZ + 0.5 * dt * k1_v;
         double posk1  = positionZ + 0.5 * dt * k1_x;
         double thetaK1 = thetaZ + 0.5 * dt * k1_theta;
         
         double k2_rho   = rho(posk1);
         double k2_x     = vk1;
         double k2_v     = -grav
                           - (0.5 / MASS) * k2_rho * drag(deploymentAngle, vk1)
                           * surfaceA(deploymentAngle)
                           * pow(vk1, 2)
                           * (1.0 / sin(thetaK1));
         double k2_theta = grav * sin(thetaK1) * cos(thetaK1) / vk1;
 
         // k3
         double vk2    = velocityZ + 0.5 * dt * k2_v;
         double posk2  = positionZ + 0.5 * dt * k2_x;
         double thetaK2 = thetaZ + 0.5 * dt * k2_theta;
         
         double k3_rho   = rho(posk2);
         double k3_x     = vk2;
         double k3_v     = -grav
                           - (0.5 / MASS) * k3_rho * drag(deploymentAngle, vk2)
                           * surfaceA(deploymentAngle)
                           * pow(vk2, 2)
                           * (1.0 / sin(thetaK2));
         double k3_theta = grav * sin(thetaK2) * cos(thetaK2) / vk2;
 
         // k4
         double vk3    = velocityZ + dt * k3_v;
         double posk3  = positionZ + dt * k3_x;
         double thetaK3 = thetaZ + dt * k3_theta;
         
         double k4_rho   = rho(posk3);
         double k4_x     = vk3;
         double k4_v     = -grav
                           - (0.5 / MASS) * k4_rho * drag(deploymentAngle, vk3)
                           * surfaceA(deploymentAngle)
                           * pow(vk3, 2)
                           * (1.0 / sin(thetaK3));
         double k4_theta = grav * sin(thetaK3) * cos(thetaK3) / vk3;
 
         // Update the state using 4th-order Runge-Kutta

         stateVector[3] += (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x) * dt / 6.0;
         stateVector[2] += (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v) * dt / 6.0;
         stateVector[6] -= (k1_theta + 2.0*k2_theta + 2.0*k3_theta + k4_theta) * dt / 6.0;

     //printf("Iterations: %d\n", iter);
     return stateVector;
}
 
