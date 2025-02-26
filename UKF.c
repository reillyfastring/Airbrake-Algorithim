#include "ukf.h"

void processModel(double sigma, double t, double* DeployAngle, double* sensorVector) {
    //drag first(1)-> velocity(1) -> position(1) -> angle(3) -> bias(3)
    for (int i = 0; i < (2*SIZE+1); i++) {
        double magA = sqrt(pow(sensorVector[1], 2) + pow(sensorVector[2], 2) + pow(sensorVector[3], 2));
        //Prediction for dragZ
        sigma[i*SIZE] += drag(DeployAngle, sigma[i*SIZE + 1]);
        //prediction for velocityZ
        sigma[i*SIZE + 1] += sin(magA - sigma[i*SIZE + 0])*t;
        //prediction for PositionZ
        sigma[i*SIZE + 2] += sigma[i*SIZE + 1]*t;
        //prediction for thetaX
        sigma[i*SIZE + 3] += sensorVector[3]*t;
        //prediction for thetaY
        sigma[i*SIZE + 4] += sensorVector[4]*t;
        //prediction for thetaZ
        sigma[i*SIZE + 5] += sensorVector[5]*t;
        //prediction for Accelerometer Bias
        sigma[i*SIZE + 6] += AccC*randN();
        //prediction for Gyroscope Bias
        sigma[i*SIZE + 7] += GyroC*randN();
        //prediction for GPS bias
        sigma[i*SIZE + 8] += GPSC*randN();
        //sensor vectors
        sigma[i*SIZE + 9] = sensorVector[0];
        sigma[i*SIZE + 10] = sensorVector[1];
        sigma[i*SIZE + 11] = sensorVector[2];
        sigma[i*SIZE + 12] = sensorVector[3];
        sigma[i*SIZE + 13] = sensorVector[4];
        sigma[i*SIZE + 14] = sensorVector[5];
        sigma[i*SIZE + 15] = sensorVector[6];
        sigma[i*SIZE + 16] = sensorVector[7];
        sigma[i*SIZE + 17] = sensorVector[8];
    }
}

void computeMean(double* sigma, double* Weights, double* stateMean) {
    for (int i = 0; i < SIZE; i++) {
        stateMean[i] = 0.0;
    }
    for (int k = 0; k < (SIZE*2+1); k++) {
        for (int j = 0; j < (SIZE); j++) {
            stateMean[j] += Weights[k]*sigma[k*SIZE + i];
        }
    }
}

void PredictCovariance(double* sigma, double* stateMean, double* PredictedCovariance, double* weights, double* ProcessNoise) {
    for (int a = 0; a < SIZE*SIZE; a++) {
        PredictedCovariance[a] = ProcessNoise[a];
    }

    for (int i = 0; i < (2*SIZE+1); i++) {
        double difference[SIZE];
        for (int k = 0; k < SIZE; k++) {
            difference[k] = sigma[i*SIZE + k] - stateMean[k];
        }

        for (int j = 0; j < SIZE; j++) {
            for (int h = 0; h < SIZE; h++) {
                PredictedCovariance[j*SIZE + h] += weights[i]*(difference[j]*difference[h]);
            }
        }
    }
}

void MeasurementFunction(double* Measurement, double* sigma) {
    for (int i = 0; i < (2*SIZE+1); i++) {
        Measurement[i*SIZE/2 + 0] = sigma[i*SIZE+9] + sigma[i*SIZE+6];
        Measurement[i*SIZE/2 + 1] = sigma[i*SIZE+10] + sigma[i*SIZE+6];
        Measurement[i*SIZE/2 + 2] = sigma[i*SIZE+11] + sigma[i*SIZE+6];
        Measurement[i*SIZE/2 + 3] = sigma[i*SIZE+12] + sigma[i*SIZE+7];
        Measurement[i*SIZE/2 + 4] = sigma[i*SIZE+13] + sigma[i*SIZE+7];
        Measurement[i*SIZE/2 + 5] = sigma[i*SIZE+14] + sigma[i*SIZE+7];
        Measurement[i*SIZE/2 + 6] = sigma[i*SIZE+15] + sigma[i*SIZE+8];
        Measurement[i*SIZE/2 + 7] = sigma[i*SIZE+16] + sigma[i*SIZE+8];
        Measurement[i*SIZE/2 + 8] = sigma[i*SIZE+17] + sigma[i*SIZE+8];
    }
}

void CrossCovariance(double* sigma, double* measurementMatrix, double* stateMean, double* measurementMean, double* crossCovariance, double* weights) {
    for (int a = 0; a < SIZE*SIZE/2; a++) {
        crossCovariance[a] = 0;
    }
    for (int b = 0; b < (2*SIZE+1); b++) {
        for (int c = 0; c < SIZE; c++) {
            double differenceState = sigma[b*SIZE + c] - stateMean[c];
            for (int d = 0; d < SIZE/2; d++) {
                double differenceMeasure = measurementMatrix[b*SIZE/2 + d] - measurementMean[d];
                crossCovariance[c*SIZE/2 + d] += weights[b]*differenceState*differenceMeasure;
            }
        }
    }
}   

void cholUpdateMulti(double* Sx, double* U, int n, int m, int sign) {
    double temp[SIZE];

    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            temp[i] = U[i*m + j];
        }
        cholUpdate(Sx, temp, n, sign);
    }
}

void cholUpdate(double* S, const double* x, int n, int sign)
{
    for (int i = 0; i < n; i++) {
        double alpha;
        if (sign > 0)  {
            alpha = Sii*Sii + xi*xi;
        }
        else if (sign < 0) {
            alpha = Sii*Sii - xi*xi;
        }

        double r = sqrt(alpha);   
        double c = r / Sii;
        double s = xi / Sii;
        S[i*n + i] = r;
        
        for (int j = i+1; j < n; j++) {
            double Sj_i = S[j*n + i];  
            double beta;
            if (sign > 0) {
                beta = Sj_i + s*x[j];
            }
            if else (sign < 0) {
                beta = Sj_i - s*x[j];
            }
            S[j*n + i] = beta/c;
            x[j] = c*x[j] - s*Sj_i;  
        }
    }
}

void updateState(double* X0, double* K, double* Zpred, double* Zinit, double* stateVector) {
    double error[SIZE/2] = {0};
    for (int i = 0; i < SIZE/2; i++) {
        error[i] = Zinit[i] - Zpred[i];
    }
    double correction[SIZE] = {0};
    matrixMultiply(K, error, correction, SIZE, SIZE/2, 1);
    for (int i = 0; i < SIZE; i++) {
        stateVector[i] = X0[i] + correction[i];
    }
}

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

void invertLowerTriangular(double* L, double* L_inv, int n) {
    for (int i = 0; i < n; i++) {
        L_inv[i*n + i] = 1.0 / L[i*n + i];
        for (int j = 0; j < i; j++) {
            double sum = 0.0;
            for (int k = j; k < i; k++) {
                sum += L[i*n + k] * L_inv[k*n + j];
            }
            L_inv[i*n + j] = -sum / L[i*n + i];
        }
        for (int j = i + 1; j < n; j++) {
            L_inv[i*n + j] = 0.0;
        }
    }
}

void multiplyTranspose(double* L_inv, double* A_inv, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += L_inv[k*n + i] * L_inv[k*n + j];
            }
            A_inv[i*n + j] = sum;
        }
    }
}

void invertMatrixCholesky(double* A, double* A_inv, int n) {
    double L[n*n];
    double L_inv[n*n];
    for (int i = 0; i < n*n; i++) {
        L[i] = A[i];
    }
    choleskyDecomp(L, n);
    invertLowerTriangular(L, L_inv, n);
    multiplyTranspose(L_inv, A_inv, n);
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

void matrixMultiply(double *A, double *B, double *C, int rowsA, int colsA, int colsB) {
    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsB; j++) {
            double sum = 0.0;
            for (int k = 0; k < colsA; k++) {
                sum += A[i * colsA + k] * B[k * colsB + j];
            }
            C[i * colsB + j] = sum;
        }
    }
}

void matrixTranspose(double* A, double* A_T, int a) {
    for (int i = 0; i < a; i++) {
        for (int j = 0; j < a; j++) {
            A_T[j * a + i] = A[i * a + j];
        }
    }
}

double randn() {
    double u1 = (rand() + 1.0) / (RAND_MAX + 1.0); 
    double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

double drag(double theta, double velocity) {
    return(0.3);
}

double rho(double height) {
    return (pow(2116*((59+459.7-.00356*(height*3.281))/518.6), 5.256)/(1718*(59+459.7-.00356*(height*3.281)))*515.4);
}

double surfaceA(angleOfDeployment) {
    surfaceA = 0.0182414692475+angleofDeployment*0.00350967*9;
}

double PredictDeploymentAngle(double* stateVector, double previousAngle) {
    double deploymentAngle = previousAngle;
    double LR = 0.001;
    double beta1 = 0.9;
    double beta2 = 0.999;
    double eps = pow(10, -8);
    int maxIterations = 1000;
    int i = 0;
    double m, v = 0;
    while(cost(deploymentAngle, stateVector) >= 25.0 && i < maxIterations) {
        double grad = computeGrad(deploymentAngle, stateVector, eps);
        m = beta1*m + (1-beta1)*grad;
        v = beta2*v + (1-beta2)*pow(grad, 2);
        double mbias = m/(1-pow(beta1, (i+1)));
        double vbias = v/(1-pow(beta2, (i+1)));
        deploymentAngle = deploymentAngle -LR*mbias/(sqrt(vbias)+eps);
        i = i + 1;
    }
}

double cost(double deploymentAngle, double* stateVector) {
    double predictedApogee = PredictApogee(deploymentAngle, stateVector);
    double distance = predictedApogee - 30000.0;
    if (distance > 0.0) {
        return pow(distance, 2);
    } else {
        return pow(distance, 2) + pow(4000.0, 2);
    }
}

double computeGrad(double deploymentAngle, double* stateVector, double epsilon) {
    double costPlus  = cost(deploymentAngle + eps, stateVector);
    double costMinus = cost(deploymentAngle - eps, stateVector);
    return (costPlus - costMinus) / (2.0 * eps);
}

double PredictApogee(double* stateVector, double deploymentAngle) {
    double velocityZ = stateVector[2];
    double positionZ = stateVector[3];
    double thetaZ = stateVector[6];
    double aZ = 0;
    double t = 0.01;
    double k1_x, k1_v, k1_theta, k1_rho, k2_x, k2_v, k2_theta, k2_rho, k3_x, k3_v, k3_theta, k3_rho, k4_x, k4_v, k4_theta, k4_rho = 0;
    int max_iterations = 0;
    while(velocityZ > 0.0 && max_iterations < 1000) {
        k1_rho = rho(positionZ);
        k1_x = velocityZ;
        k1_v = -grav -(0.5/MASS)*k1_rho*drag(velocityZ, deploymentAngle)*surfaceA(deploymentAngle)*pow(velocityZ, 2)*(1.0/math.sin(thetaZ)));
        k1_theta = grav*math.sin(thetaZ)*math.cos(thetaZ)/velocityZ;

        k2_rho = rho(positionZ + velocityZ*t/2.0);
        k2_x = velocityZ+k1_v*t/2.0;
        k2_v = -grav -(0.5/MASS)*k2_rho*drag((velocityZ+k1_v*t/2.0), deploymentAngle)*surfaceA(deploymentAngle)*pow((velocityZ+k1_v*t/2.0), 2)*(1.0/math.sin(thetaZ+k1_theta*t/2.0)));
        k2_theta = grav*math.sin(thetaZ+k1_theta*t/2.0)*math.cos(thetaZ+k1_theta*t/2.0)/(velocityZ+k1_v*t/2.0);

        k3_rho = rho(positionZ + velocityZ*t/2);
        k3_x = positionZ*(velocityZ+k2_v*t/2);
        k3_v = -grav -(1/2*MASS)*k2_rho*drag((velocityZ+k2_v*t/2), deploymentAngle)*pow((velocityZ+k2_v*t/2), 2)*(1/math.sin(thetaZ+k2_theta*t/2)));

        k3_rho = rho(positionZ + velocityZ*t/2.0);
        k3_x = velocityZ+k2_v*t/2.0;
        k3_v = -grav -(0.5/MASS)*k3_rho*drag((velocityZ+k2_v*t/2.0), deploymentAngle)*surfaceA(deploymentAngle)*pow((velocityZ+k2_v*t/2.0), 2)*(1.0/math.sin(thetaZ+k2_theta*t/2.0)));
        k3_theta = grav*math.sin(thetaZ+k2_theta*t/2)*math.cos(thetaZ+k2_theta*t/2)/(velocityZ+k2_v*t/2);

        k4_rho = rho(positionZ + velocityZ*t);
        k4_x = velocityZ+k3_v*t;
        k4_v = -grav -(0.5/MASS)*k4_rho*drag((velocityZ+k3_v*t), deploymentAngle)*surfaceA(deploymentAngle)*pow((velocityZ+k3_v*t), 2)*(1.0/math.sin(thetaZ+k3_theta*t)));
        k4_theta = grav*math.sin(thetaZ+k3_theta*t)*math.cos(thetaZ+k3_theta*t)/(velocityZ+k3_v*t);
        
        positionZ = positionZ + (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x)*t/6.0;
        velocityZ = velocityZ + (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v)*t/6.0;
        thetaZ = thetaZ + (k1_theta + 2.0*k2_theta + 2.0*k3_theta + k4_theta)*t/6.0;
        max_iterations = max_iterations + 1;
    }
    return(positionZ)
}