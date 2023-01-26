#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

// sigmoid function
double sigmoid(double x) {
    return 1 / (1 + exp(-x));
}

// derivative of sigmoid function
double sigmoidDerivative(double x) {
    return sigmoid(x) * (1 - sigmoid(x));
}

int main() {
    // neural network configuration
    static const int numInputs = 2;
    static const int numHidden = 2;
    static const int numOutputs = 2;

    double hiddenLayer[numHidden];
    double outputLayer[numOutputs];

    // weights
    double weightsInputToHidden[numInputs][numHidden];
    double weightsHiddenToOutputs[numHidden][numOutputs];

    // bias
    double biasHidden[numHidden];
    double biasOutput[numOutputs];

    // training data
    const int totalTrainingData = 4;

    double trainingInputs[totalTrainingData][numInputs] = {
        {0, 0},
        {0, 1},
        {1, 0},
        {1, 1}
    };

    double correctTrainingOutputs[totalTrainingData][numOutputs] = {
        {1, 0},
        {0, 1},
        {0, 1},
        {1, 0}
    };

    // initialize all values to zero instead of randomizing so the output is the same every time
    // initialize hidden weights and bias to defaultWeights

    double defaultWeight1 = 0.420;
    double defaultWeight2 = 0.69;

    for (int i = 0; i < numHidden; i++) {

        // initialize to weights 1 or 2 i if i is even or odd
        biasHidden[i] = i % 2 == 0 ? defaultWeight1 : defaultWeight2;

        for (int j = 0; j < numInputs; j++) {
            weightsInputToHidden[j][i] = (j + i) % 2 == 0 ? defaultWeight1 : defaultWeight2;
        }
    }

    for (int i = 0; i < numOutputs; i++) {

        // initialize to weights 1 or 2 i if i is even or odd
        biasOutput[i] = i % 2 == 0 ? defaultWeight1 : defaultWeight2;
        for (int j = 0; j < numHidden; j++) {
            weightsHiddenToOutputs[j][i] = (j + i) % 2 == 0 ? defaultWeight1 : defaultWeight2;
        }
    }

    // initialize hidden weights and bias to random values if you want...

    // for (int i = 0; i < numHidden; i++) {

    //     // initialize to random values
    //     biasHidden[i] = (double)rand() / RAND_MAX;

    //     for (int j = 0; j < numInputs; j++) {
    //         weightsInputToHidden[j][i] = (double)rand() / RAND_MAX;
    //     }
    // }

    // for (int i = 0; i < numOutputs; i++) {

    //     // initialize to random values
    //     biasOutput[i] = (double)rand() / RAND_MAX;
    //     for (int j = 0; j < numHidden; j++) {
    //         weightsHiddenToOutputs[j][i] = (double)rand() / RAND_MAX;
    //     }
    // }

    printf("\n");
    printf("STEP_0: Initialize weights and bias:\n");

    // print weights
    printf("---------------\n");
    printf("weightsInputToHidden: \n");
    for (int i = 0; i < numInputs; i++) {
        for (int j = 0; j < numHidden; j++) {
            printf("\t [%d][%d] %f ", i, j, weightsInputToHidden[i][j]);
        }
        printf("\n");
    }

    // print bias
    printf("---------------\n");
    printf("biasHidden: \n");
    for (int i = 0; i < numHidden; i++) {
        printf("\t%f ", biasHidden[i]);
    }
    printf("\n");
    printf("\n");

    printf("---------------\n");
    printf("weightsHiddenToOutputs: \n");
    for (int i = 0; i < numHidden; i++) {
        for (int j = 0; j < numOutputs; j++) {
            printf("\t [%d][%d] %f ", i, j, weightsHiddenToOutputs[i][j]);
        }
        printf("\n");
    }

    // print bias
    printf("---------------\n");
    printf("biasOutput: \n");
    for (int i = 0; i < numOutputs; i++) {
        printf("\t%f ", biasOutput[i]);
    }
    printf("\n");

    // simple forward propagation before training
    printf("\n");
    printf("STEP_1: Calculate accuracy before training:\n");
    printf("---------------\n");
    int correctMatchedBeforeTraining = 0;
    for (int i = 0; i < totalTrainingData; i++) {

        // print input data
        printf("training data: %d\n", i);
        for (int j = 0; j < numInputs; j++) {
            printf("\t%d: %f \n", j, trainingInputs[i][j]);
        }
        printf("\n");

        // calculate hidden layer
        double computedHiddenLayer[numHidden];

        for (int j = 0; j < numHidden; j++) {
            double sum = 0;
            for (int k = 0; k < numInputs; k++) {
                sum += trainingInputs[i][k] * weightsInputToHidden[k][j];
            }
            sum += biasHidden[j];
            computedHiddenLayer[j] = sigmoid(sum);
        }

        // calculate output layer
        double computedOutputLayer[numOutputs];

        for (int j = 0; j < numOutputs; j++) {
            double sum = 0;
            for (int k = 0; k < numHidden; k++) {
                sum += computedHiddenLayer[k] * weightsHiddenToOutputs[k][j];
            }
            sum += biasOutput[j];
            computedOutputLayer[j] = sigmoid(sum);
        }

        // find max value in computedOutputLayer
        double max = 0;
        int maxIndex = 0;
        for (int j = 0; j < numOutputs; j++) {
            if (computedOutputLayer[j] > max) {
                max = computedOutputLayer[j];
                maxIndex = j;
            }
        }

        double predictedOutput[numOutputs];

        // set all values to zero
        for (int j = 0; j < numOutputs; j++) {
            predictedOutput[j] = 0;

            // if maxIndex is equal to j, set predictedOutput[j] to 1
            if (maxIndex == j) {
                predictedOutput[j] = 1;
            }
        }

        // print predicted output
        printf("\tpredicted output: \n");
        for (int j = 0; j < numOutputs; j++) {
            printf("\t\t%d: %f \n", j, predictedOutput[j]);
        }

        double correctOutput[numOutputs];

        // get correct output from correctTrainingOutputs
        for (int j = 0; j < numOutputs; j++) {
            correctOutput[j] = correctTrainingOutputs[i][j];
        }
        
        bool correct = false;

        // check if predicted output is equal to correct output
        for (int j = 0; j < numOutputs; j++) {
            if (predictedOutput[j] == correctOutput[j]) {
                correct = true;
            } else {
                correct = false;
                break;
            }
        }

        // print if predicted output is correct
        if (correct) {
            printf("\tcorrect\n");
            correctMatchedBeforeTraining++;
        } else {
            printf("\tincorrect\n");
        }

        printf("\n");
    }

    // print accuracy
    printf("before training accuracy: %d%%\n", (int)((double) correctMatchedBeforeTraining * 100 / totalTrainingData));
    printf("\n");

    // training step
    int epoch = 7;
    double learningRate = 0.1;

    printf("STEP_2: Begin training...\n");
    printf("---------------\n");

    for (int i = 0; i < epoch; i++) {

        double totalError = 0;
        for (int j = 0; j < totalTrainingData; j++) {

            printf("\n");
            printf("training data: %d\n", j);
            printf("\n");
            // print input data
            printf("input data: \n");
            for (int k = 0; k < numInputs; k++) {
                printf("\t%d: %f \n", k, trainingInputs[j][k]);
            }
            printf("\n");

            // do forward propagation

            // calculate hidden layer
            double computedHiddenLayer[numHidden];

            for (int k = 0; k < numHidden; k++) {
                double sum = 0;
                for (int l = 0; l < numInputs; l++) {
                    sum += trainingInputs[j][l] * weightsInputToHidden[l][k];
                }
                sum += biasHidden[k];
                computedHiddenLayer[k] = sigmoid(sum);
            }

            // calculate output layer
            double computedOutputLayer[numOutputs];

            for (int k = 0; k < numOutputs; k++) {
                double sum = 0;
                for (int l = 0; l < numHidden; l++) {
                    sum += computedHiddenLayer[l] * weightsHiddenToOutputs[l][k];
                }
                sum += biasOutput[k];
                computedOutputLayer[k] = sigmoid(sum);
            }

            // ok! Now it's time for back propagation

            for (int k = 0; k < numOutputs; k++) {
                double errorRate = pow((correctTrainingOutputs[j][k] - computedOutputLayer[k]), 3);

                printf("\t\toutputNum: %d\n", k);

                printf("\t\tcorrectTrainingOutputs (outputNum: %d): %f\n", k, correctTrainingOutputs[j][k]);
                printf("\t\tcomputedOutputLayer (outputNum: %d): %f\n", k, computedOutputLayer[k]);
                printf("\t\terrorRate %f\n\n", errorRate);

                totalError += fabs(errorRate);

                // get partial derivative for this output with represent to each hidden input
                for (int l = 0; l < numHidden; l++) {
                    printf("\t\t\thiddenNum: %d\n", l);

                    printf("\t\t\tweightsHiddenToOutputs (hiddenNum: %d, outputNum: %d): %f\n", l, k, weightsHiddenToOutputs[l][k]);

                    double weightsHiddenToOutputsPartialDerivative = sigmoidDerivative(weightsHiddenToOutputs[l][k]);
                    printf("\t\t\tweightsHiddenToOutputsPartialDerivative (hiddenNum: %d, outputNum: %d): %f\n", l, k, weightsHiddenToOutputsPartialDerivative);
                
                    double weightsHiddenToOutputsPartialDerivativeWithErrorMagnitude = weightsHiddenToOutputsPartialDerivative * errorRate;
                    printf("\t\t\tweightsHiddenToOutputsPartialDerivativeWithErrorMagnitude (hiddenNum: %d, outputNum: %d): %f\n", l, k, weightsHiddenToOutputsPartialDerivativeWithErrorMagnitude);
                    printf("\n");

                    // get partial derivative for this hidden output with represent to each input
                    for(int m = 0; m < numInputs; m++) {
                        printf("\t\t\t\tinputNum: %d\n", m);

                        printf("\t\t\t\tweightsInputToHidden (inputNum: %d, hiddenNum: %d): %f\n", m, l, weightsInputToHidden[m][l]);

                        double weightsInputToHiddenPartialDerivative = sigmoidDerivative(weightsInputToHidden[m][l]);
                        printf("\t\t\t\tweightsInputToHiddenPartialDerivative (inputNum: %d, hiddenNum: %d): %f\n", m, l, weightsInputToHiddenPartialDerivative);

                        double weightsInputToHiddenPartialDerivativeWithErrorMagnitude = weightsInputToHiddenPartialDerivative * weightsHiddenToOutputsPartialDerivativeWithErrorMagnitude;
                        printf("\t\t\t\tweightsInputToHiddenPartialDerivativeWithErrorMagnitude (inputNum: %d, hiddenNum: %d): %f\n", m, l, weightsInputToHiddenPartialDerivativeWithErrorMagnitude);
                        printf("\n");

                        // update weights input to hidden to reflect new derivative
                        weightsInputToHidden[m][l] += weightsInputToHiddenPartialDerivativeWithErrorMagnitude;

                        // get partial derivative for this hidden output with represent to each bias
                        double biasHiddenPartialDerivative = sigmoidDerivative(biasHidden[l]);
                        printf("\t\t\t\tbiasHiddenPartialDerivative (hiddenNum: %d): %f\n", l, biasHiddenPartialDerivative);

                        double biasHiddenPartialDerivativeWithErrorMagnitude = biasHiddenPartialDerivative * weightsHiddenToOutputsPartialDerivativeWithErrorMagnitude;
                        printf("\t\t\t\tbiasHiddenPartialDerivativeWithErrorMagnitude (hiddenNum: %d): %f\n", l, biasHiddenPartialDerivativeWithErrorMagnitude);
                        printf("\n");

                        // update bias hidden to reflect new derivative
                        biasHidden[l] += biasHiddenPartialDerivativeWithErrorMagnitude;
                    }


                    // update weights hidden to output to reflect new derivative
                    weightsHiddenToOutputs[l][k] += weightsHiddenToOutputsPartialDerivativeWithErrorMagnitude;

                    // get partial derivative for this output with represent to each bias
                    double biasOutputPartialDerivative = sigmoidDerivative(biasOutput[k]);
                    printf("\t\t\tbiasOutputPartialDerivative (outputNum: %d): %f\n", k, biasOutputPartialDerivative);

                    double biasOutputPartialDerivativeWithErrorMagnitude = biasOutputPartialDerivative * errorRate;
                    printf("\t\t\tbiasOutputPartialDerivativeWithErrorMagnitude (outputNum: %d): %f\n", k, biasOutputPartialDerivativeWithErrorMagnitude);
                    printf("\n");

                    // update bias output to reflect new derivative
                    biasOutput[k] += biasOutputPartialDerivativeWithErrorMagnitude;
                }

                printf("\n");
            }

        }

        // print weights
        printf("---------------\n");
        printf("weightsInputToHidden: \n");
        for (int i = 0; i < numInputs; i++) {
            for (int j = 0; j < numHidden; j++) {
                printf("\t [%d][%d] %f ", i, j, weightsInputToHidden[i][j]);
            }
            printf("\n");
        }

        // print bias
        printf("---------------\n");
        printf("biasHidden: \n");
        for (int i = 0; i < numHidden; i++) {
            printf("\t%f ", biasHidden[i]);
        }
        printf("\n");
        printf("\n");

        printf("---------------\n");
        printf("weightsHiddenToOutputs: \n");
        for (int i = 0; i < numHidden; i++) {
            for (int j = 0; j < numOutputs; j++) {
                printf("\t [%d][%d] %f ", i, j, weightsHiddenToOutputs[i][j]);
            }
            printf("\n");
        }

        // print bias
        printf("---------------\n");
        printf("biasOutput: \n");
        for (int i = 0; i < numOutputs; i++) {
            printf("\t%f ", biasOutput[i]);
        }
        printf("\n");

        printf("epoch: %d errorRate: %f\n", i, totalError);
    }
    printf("\n");
    printf("TRAINING COMPLETE! If all goes well, you should see the error rate gradually decreasing to zero as the epoch increases.\n\n");

    // simple forward propagation after training
    printf("\n");
    printf("STEP_3: Calculate accuracy after training:\n");
    printf("---------------\n");
    int correctMatchedAfterTraining = 0;
    for (int i = 0; i < totalTrainingData; i++) {

        // print input data
        printf("training data: %d\n", i);
        for (int j = 0; j < numInputs; j++) {
            printf("\t%d: %f \n", j, trainingInputs[i][j]);
        }
        printf("\n");

        // calculate hidden layer
        double computedHiddenLayer[numHidden];

        for (int j = 0; j < numHidden; j++) {
            double sum = 0;
            for (int k = 0; k < numInputs; k++) {
                sum += trainingInputs[i][k] * weightsInputToHidden[k][j];
            }
            sum += biasHidden[j];
            computedHiddenLayer[j] = sigmoid(sum);
        }

        // calculate output layer
        double computedOutputLayer[numOutputs];

        for (int j = 0; j < numOutputs; j++) {
            double sum = 0;
            for (int k = 0; k < numHidden; k++) {
                sum += computedHiddenLayer[k] * weightsHiddenToOutputs[k][j];
            }
            sum += biasOutput[j];
            computedOutputLayer[j] = sigmoid(sum);
        }

        // find max value in computedOutputLayer
        double max = 0;
        int maxIndex = 0;
        for (int j = 0; j < numOutputs; j++) {
            if (computedOutputLayer[j] > max) {
                max = computedOutputLayer[j];
                maxIndex = j;
            }
        }

        double predictedOutput[numOutputs];

        // set all values to zero
        for (int j = 0; j < numOutputs; j++) {
            predictedOutput[j] = 0;

            // if maxIndex is equal to j, set predictedOutput[j] to 1
            if (maxIndex == j) {
                predictedOutput[j] = 1;
            }
        }

        // print predicted output
        printf("\tpredicted output: \n");
        for (int j = 0; j < numOutputs; j++) {
            printf("\t\t%d: %f \n", j, predictedOutput[j]);
        }

        double correctOutput[numOutputs];

        // get correct output from correctTrainingOutputs
        for (int j = 0; j < numOutputs; j++) {
            correctOutput[j] = correctTrainingOutputs[i][j];
        }
        
        bool correct = false;

        // check if predicted output is equal to correct output
        for (int j = 0; j < numOutputs; j++) {
            if (predictedOutput[j] == correctOutput[j]) {
                correct = true;
            } else {
                correct = false;
                break;
            }
        }

        // print if predicted output is correct
        if (correct) {
            printf("\tcorrect\n");
            correctMatchedAfterTraining++;
        } else {
            printf("\tincorrect\n");
        }

        printf("\n");
    }

    // print accuracy
    printf("after training accuracy: %d%%\n", (int)((double) correctMatchedAfterTraining * 100 / totalTrainingData));
    printf("\n");

}