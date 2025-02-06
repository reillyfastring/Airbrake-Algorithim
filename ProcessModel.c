#define MASS 20
#define grav 9.80665
#define AccC 0.01
#define GyroC 0.001
#define GPSC 0.05

void processModel(double* A, int n, double* t, double* DeployAngle, double* sensorVector) {
    //drag first(1)-> velocity(1) -> position(1) -> angle(3) -> bias(3)
    for (int i = 1; i < (n+1); i++) {
        magA = sqrt(pow(sensorVector[1], 2) + pow(sensorVector[2], 2) + pow(sensorVector[3], 2));
        //Prediction for dragZ
        sigma[i*SIZE + 1] += drag(theta, sigma[i*SIZE + 2]);
        //prediction for velocityZ
        sigma[i*SIZE + 2] += sin(magA - sigma[i*SIZE + 1])*t;
        //prediction for PositionZ
        sigma[i*SIZE + 3] += sigma[i*SIZE + 2]*t;
        //prediction for thetaX
        sigma[i*SIZE + 4] += sensorVector[4]*t;
        //prediction for thetaY
        sigma[i*SIZE + 5] += sensorVector[5]*t;
        //prediction for thetaZ
        sigma[i*SIZE + 6] += sensorVector[6]*t;
        //prediction for Accelerometer Bias
        sigma[i*SIZE + 7] += AccC*randN();
        //prediction for Gyroscope Bias
        sigma[i*SIZE + 8] += GryoC*randN();
        //prediction for GPS bias
        sigma[i*SIZE + 9] += GPSC*randN();
    }
}

