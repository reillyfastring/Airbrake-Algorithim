void computeMean(double* sigma, double* Weights, double* stateMean) {
    for (int i = 0; i < SIZE; i++) {
        for (int k = 0; k < (SIZE*2+1); k++) {
            stateMean[i] += Weights[2*k+2]*sigma[i*SIZE + k];
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
            difference[k] = sigma[i*SIZE + i] - stateMean[i];
        }

        for (int j = 0; j < SIZE; j++) {
            for (int h = 0; h < SIZE; h++) {
                PredictedCovariance[j*SIZE + h] += weights[i]*(difference[j]*difference[h]);
            }
        }
    }
}

void MeasurementFunction(double* stateVector, double* Measrurement) {
    
}

