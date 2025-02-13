#include "UKF.c"
#include "PredictApogee.c"

double PredictDeploymentAngle(double stateVector, double previousAngle) {
    double distance = 0;
    double cost = 0;
    while(cost >= 25, iterations >= 100) {
        distance = (PredictApogee(previousAngle, stateVector) - 30000);
        if (distance > 0) {
            cost = pow(distance, 2);
        }
        else {
            cost = pow(distance, 2) + pow(4000, 2);
        }
        
    }
}

double PredictApogee(double* stateVector, double deploymentAngle) {
    double velocityZ = stateVector[2];
    double positionZ = stateVector[3];
    double thetaZ = stateVector[6];
    double aZ = 0;
    double t = 0.01;
    double k1_x, k1_v, k1_theta, k1_rho, k2_x, k2_v, k2_theta, k2_rho, k3_x, k3_v, k3_theta, k3_rho, k4_x, k4_v, k4_theta, k4_rho = 0;
    int max_iterations = 0;
    while(velocityZ > 0 || max_iterations > 1000) {
        k1_rho = rho(positionZ);
        k1_x = positionZ*velocityZ;
        k1_v = -grav -(1/2*MASS)*k1_rho*drag(velocityZ, deploymentAngle)*pow(velocityZ, 2)*(1/math.sin(thetaZ)));
        k1_theta = grav*math.sin(thetaZ)*math.cos(thetaZ)/velocityZ;

        k2_rho = rho(positionZ + velocityZ*t/2);
        k2_x = positionZ*(velocityZ+k1_v*t/2);
        k2_v = -grav -(1/2*MASS)*k1_rho*drag((velocityZ+k1_v*t/2), deploymentAngle)*pow((velocityZ+k1_v*t/2), 2)*(1/math.sin(thetaZ+k1_theta*t/2)));
        k2_theta = grav*math.sin(thetaZ+k1_theta*t/2)*math.cos(thetaZ+k1_theta*t/2)/(velocityZ+k1_v*t/2);


        k3_rho = rho(positionZ + velocityZ*t/2);
        k3_x = positionZ*(velocityZ+k2_v*t/2);
        k3_v = -grav -(1/2*MASS)*k2_rho*drag((velocityZ+k2_v*t/2), deploymentAngle)*pow((velocityZ+k2_v*t/2), 2)*(1/math.sin(thetaZ+k2_theta*t/2)));
        k3_theta = grav*math.sin(thetaZ+k2_theta*t/2)*math.cos(thetaZ+k2_theta*t/2)/(velocityZ+k2_v*t/2);

        k4_rho = rho(positionZ + velocityZ*t);
        k4_x = positionZ*(velocityZ+k3_v*t);
        k4_v = -grav -(1/2*MASS)*k3_rho*drag((velocityZ+k3_v*t), deploymentAngle)*pow((velocityZ+k3_v*t), 2)*(1/math.sin(thetaZ+k3_theta*t)));
        k4_theta = grav*math.sin(thetaZ+k3_theta*t)*math.cos(thetaZ+k3_theta*t)/(velocityZ+k3_v*t);
        
        positionZ = positionZ + (k1_x + 2*k2_x + 2*k3_x + k4_x)*t/6;
        velocityZ = velocityZ + (k1_v + 2*k2_v + 2*k3_v + k4_v)*t/6;
        thetaZ = thetaZ + (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta)*t/6;
        max_iterations = max_iterations + 1;
    }
    return(positionZ)
}

