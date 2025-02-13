#include "UKF.c"
#include "PredictApogee.c"

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

