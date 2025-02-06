#define MASS 20
#define grav

void processModel(double* A, int n, double* t, double* DeployAngle) {
    //drag first(1)-> velocity(3) -> position(3) -> angle(3) -> bias(3)
    for (int i = 1; i < (n+1); i++) {
        double magV = sqrt(pow(sigma[i*SIZE + 2], 2) + pow(sigma[i*SIZE + 3], 2) + pow(sigma[i*SIZE + 4], 2))
        sigma[i*SIZE + 1] += t*drag(magV, DeployAngle);
        k1 = -grav - 1/(2*MASS)*
        sigma[i*SIZE + 2] += 
    }
}