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

// given an error rate function of (correct - predicted) ^ 2
// this is the partial derivative of the error rate wrt to the 'x' value inside a sigmoid function passed alongside the predicted value
// Given this is a partial derivative: 
//      a is the correct value
//      b is the summed outputs of the other neurons
//      c is the direct input
double errorRateDerivativeWRTXInsideSigmoid(double x, double a, double b, double c) {
    return ( (2 * c * exp(-c * x * -b )) * ( ( a * (1 + exp(-b - (c * x) ) ) ) - 1) ) / pow((1 + exp(-b - (c * x)) ), 3);
}

int main() {
    // neural network configuration
    static const int numInputs = 2;
    static const int numHidden = 2;
    static const int numOutputs = 2;

    // weights
    double hiddenWeightsForEachInputNeuron[numHidden][numInputs];
    double outputWeightsForEachHiddenNeuron[numOutputs][numHidden];

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
            hiddenWeightsForEachInputNeuron[j][i] = (j + i) % 2 == 0 ? defaultWeight1 : defaultWeight2;
        }
    }

    for (int i = 0; i < numOutputs; i++) {

        // initialize to weights 1 or 2 i if i is even or odd
        biasOutput[i] = i % 2 == 0 ? defaultWeight1 : defaultWeight2;
        for (int j = 0; j < numHidden; j++) {
            outputWeightsForEachHiddenNeuron[j][i] = (j + i) % 2 == 0 ? defaultWeight1 : defaultWeight2;
        }
    }

    // initialize hidden weights and bias to random values if you want...

    // for (int i = 0; i < numHidden; i++) {

    //     // initialize to random values
    //     biasHidden[i] = (double)rand() / RAND_MAX;

    //     for (int j = 0; j < numInputs; j++) {
    //         hiddenWeightsForEachInputNeuron[j][i] = (double)rand() / RAND_MAX;
    //     }
    // }

    // for (int i = 0; i < numOutputs; i++) {

    //     // initialize to random values
    //     biasOutput[i] = (double)rand() / RAND_MAX;
    //     for (int j = 0; j < numHidden; j++) {
    //         outputWeightsForEachHiddenNeuron[j][i] = (double)rand() / RAND_MAX;
    //     }
    // }

    printf("\n");
    printf("STEP_0: Initialize weights and bias:\n");

    // print weights
    printf("---------------\n");
    printf("hiddenWeightsForEachInputNeuron: \n");
    for (int i = 0; i < numInputs; i++) {
        for (int j = 0; j < numHidden; j++) {
            printf("\t [%d][%d] %f ", i, j, hiddenWeightsForEachInputNeuron[i][j]);
        }
        printf("\n");
    }

    // print bias
    printf("---------------\n");
    printf("biasHidden: \n");
    for (int i = 0; i < numHidden; i++) {
        printf("\t biasHidden[%d] %f ", i, biasHidden[i]);
    }
    printf("\n");
    printf("\n");

    printf("---------------\n");
    printf("outputWeightsForEachHiddenNeuron: \n");
    for (int i = 0; i < numHidden; i++) {
        for (int j = 0; j < numOutputs; j++) {
            printf("\t [%d][%d] %f ", i, j, outputWeightsForEachHiddenNeuron[i][j]);
        }
        printf("\n");
    }

    // print bias
    printf("---------------\n");
    printf("biasOutput: \n");
    for (int i = 0; i < numOutputs; i++) {
        printf("\t biasOutput[%d] %f ", i, biasOutput[i]);
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

        // calculate hidden outputs
        double computedHiddenNeuronOutput[numHidden];

        for (int j = 0; j < numHidden; j++) {
            double sum = 0;
            for (int k = 0; k < numInputs; k++) {
                sum += trainingInputs[i][k] * hiddenWeightsForEachInputNeuron[j][k];
            }
            sum += biasHidden[j];
            computedHiddenNeuronOutput[j] = sigmoid(sum);
        }

        // calculate output outputs
        double computedOutputNeuronOutput[numOutputs];

        for (int j = 0; j < numOutputs; j++) {
            double sum = 0;
            for (int k = 0; k < numHidden; k++) {
                sum += computedHiddenNeuronOutput[k] * outputWeightsForEachHiddenNeuron[j][k];
            }
            sum += biasOutput[j];
            computedOutputNeuronOutput[j] = sigmoid(sum);
        }

        // find max value in computedOutputNeuronOutput
        double max = 0;
        int maxIndex = 0;
        for (int j = 0; j < numOutputs; j++) {
            if (computedOutputNeuronOutput[j] > max) {
                max = computedOutputNeuronOutput[j];
                maxIndex = j;
            }
        }

        // print computed output
        printf("\tcomputed output: \n");
        for (int j = 0; j < numOutputs; j++) {
            printf("\t\t%d: %f \n", j, computedOutputNeuronOutput[j]);
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
    int epoch = 1000;
    double learningRate = 0.01;

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

            // calculate hidden neuron output values
            double computedHiddenNeuronOutput[numHidden];

            for (int k = 0; k < numHidden; k++) {
                double sum = 0;
                for (int l = 0; l < numInputs; l++) {
                    sum += trainingInputs[j][l] * hiddenWeightsForEachInputNeuron[k][l];
                }
                sum += biasHidden[k];
                computedHiddenNeuronOutput[k] = sigmoid(sum);
            }

            // calculate output neuron output values
            double computedOutputNeuronOutput[numOutputs];

            for (int k = 0; k < numOutputs; k++) {
                double sum = 0;
                for (int l = 0; l < numHidden; l++) {
                    sum += computedHiddenNeuronOutput[l] * outputWeightsForEachHiddenNeuron[k][l];
                }
                sum += biasOutput[k];
                computedOutputNeuronOutput[k] = sigmoid(sum);
            }

            // ok! Now it's time for back propagation

            double outputWeightsForEachHiddenNeuronGradient[numOutputs][numHidden];
            double biasOutputGradient[numOutputs];
            double hiddenWeightsForEachInputNeuronGradient[numHidden][numInputs];
            double biasHiddenGradient[numHidden];

            // initialize all gradients to zero
            for (int k = 0; k < numOutputs; k++) {
                for (int l = 0; l < numHidden; l++) {
                    outputWeightsForEachHiddenNeuronGradient[k][l] = 0;
                }
                biasOutputGradient[k] = 0;
            }

            for (int k = 0; k < numHidden; k++) {
                for (int l = 0; l < numInputs; l++) {
                    hiddenWeightsForEachInputNeuronGradient[k][l] = 0;
                }
                biasHiddenGradient[k] = 0;
            }

            // calculate error for each output neuron and update the gradient
            for (int k = 0; k < numOutputs; k++) {

                double deltaPredictedFromCorrect = (correctTrainingOutputs[j][k] - computedOutputNeuronOutput[k]);

                // the cost/loss function as the cool kids call it... I prefer error rate
                double errorRate = pow(deltaPredictedFromCorrect, 2);

                printf("\n");
                printf("\t\toutputNum: %d\n", k);

                printf("\t\tcorrectTrainingOutputs (outputNum: %d): %f\n", k, correctTrainingOutputs[j][k]);
                printf("\t\tcomputedOutputNeuronOutput (outputNum: %d): %f\n", k, computedOutputNeuronOutput[k]);
                printf("\t\terrorRate %f\n", errorRate);

                totalError += errorRate;

                // calculate partial derivative of the error rate with respect to the weight of each hidden neuron
                for (int l = 0; l < numHidden; l++) {
                    printf("\n");
                    printf("\t\t\thiddenNum: %d\n", l);
                    printf("\t\t\toutputWeightsForThisHiddenNeuron (outputNum: %d, hiddenNum: %d): %f\n", k, l, outputWeightsForEachHiddenNeuron[k][l]);
                    printf("\n");
                    double xParameter = outputWeightsForEachHiddenNeuron[k][l];
                    double aParameter = correctTrainingOutputs[j][k];
                    double bParameter = 0.0;
                    double cParameter = computedHiddenNeuronOutput[l];

                    for (int ll = 0; ll < numHidden; ll++) {
                        if(ll != l){
                            bParameter += computedHiddenNeuronOutput[ll] * outputWeightsForEachHiddenNeuron[k][ll];
                        }
                    }
                    bParameter += biasOutput[k];
                    

                    double partialDerivativeErrorRateWrtHiddenNeuronWeight = learningRate * errorRateDerivativeWRTXInsideSigmoid(xParameter, aParameter, bParameter, cParameter);
                    printf("\t\t\tpartialDerivativeErrorRateWrtHiddenNeuronWeight %f \n", partialDerivativeErrorRateWrtHiddenNeuronWeight);
                    printf("\t\t\t\t parameter a: %f parameter b: %f parameter c: %f\n", aParameter, bParameter, cParameter);
                    // update outputWeightsForEachHiddenNeuron with gradient
                    outputWeightsForEachHiddenNeuronGradient[k][l] += partialDerivativeErrorRateWrtHiddenNeuronWeight;

                    // calculate partial derivative of the error rate with respect to the weight of each input neuron
                    for (int m = 0; m < numInputs; m++) {
                        printf("\n");
                        printf("\t\t\t\t\tinputNum: %d\n", m);
                        printf("\t\t\t\t\thiddenWeightsForThisInputNeuron (hiddenNum: %d, inputNum: %d): %f\n", l, m, hiddenWeightsForEachInputNeuron[l][m]);
                        printf("\n");
                        double xParameter = hiddenWeightsForEachInputNeuron[l][m];
                        double aParameter = correctTrainingOutputs[j][k];
                        double bParameter = 0.0;
                        double cParameter = trainingInputs[j][m];

                        for (int mm = 0; mm < numInputs; mm++) {
                            if(mm != m){
                                bParameter += trainingInputs[j][mm] * hiddenWeightsForEachInputNeuron[l][mm];
                            }
                        }

                        double partialDerivativeErrorRateWrtInputNeuronWeight = learningRate * errorRateDerivativeWRTXInsideSigmoid(xParameter, aParameter, bParameter, cParameter);
                        printf("\t\t\t\t\tpartialDerivativeErrorRateWrtInputNeuronWeight %f \n", partialDerivativeErrorRateWrtInputNeuronWeight);
                        printf("\t\t\t\t\t\t parameter a: %f parameter b: %f parameter c: %f\n", aParameter, bParameter, cParameter);
                        // update hiddenWeightsForEachInputNeuron with gradient
                        hiddenWeightsForEachInputNeuronGradient[l][m] += partialDerivativeErrorRateWrtInputNeuronWeight;
                    }

                    printf("\n");

                    // calculate partial derivative of the error rate with respect to the bias of the hidden neuron
                    printf("\t\t\t\tbiasHidden (hiddenNum: %d): %f\n", l, biasHidden[l]);

                    double xParameterBias = biasHidden[l];
                    double aParameterBias = correctTrainingOutputs[j][k];
                    double bParameterBias = 0.0;
                    double cParameterBias = 1.0;

                    for (int m = 0; m < numInputs; m++) {
                        bParameterBias += trainingInputs[j][m] * hiddenWeightsForEachInputNeuron[l][m];
                    }

                    double partialDerivativeErrorRateWrtBiasHidden = learningRate * errorRateDerivativeWRTXInsideSigmoid(xParameterBias, aParameterBias, bParameterBias, cParameterBias);
                    printf("\t\t\t\tpartialDerivativeErrorRateWrtBiasHidden %f \n", partialDerivativeErrorRateWrtBiasHidden);
                    printf("\t\t\t\t\t parameter a: %f parameter b: %f parameter c: %f\n", aParameterBias, bParameterBias, cParameterBias);
                    // update biasHidden with gradient
                    biasHiddenGradient[l] += partialDerivativeErrorRateWrtBiasHidden;
                }

                printf("\n");

                // calculate partial derivative of the error rate with respect to the bias of the output neuron
                printf("\t\tbiasOutput (outputNum: %d): %f\n", k, biasOutput[k]);

                double xParameter = biasOutput[k];
                double aParameter = correctTrainingOutputs[j][k];
                double bParameter = 0.0;
                double cParameter = 1.0;

                for (int l = 0; l < numHidden; l++) {
                    bParameter += computedHiddenNeuronOutput[l] * outputWeightsForEachHiddenNeuron[k][l];
                }

                double partialDerivativeErrorRateWrtBiasOutput = learningRate * errorRateDerivativeWRTXInsideSigmoid(xParameter, aParameter, bParameter, cParameter);
                printf("\t\tpartialDerivativeErrorRateWrtBiasOutput %f \n", partialDerivativeErrorRateWrtBiasOutput);
                printf("\t\t\t parameter a: %f parameter b: %f parameter c: %f\n", aParameter, bParameter, cParameter);
                // update biasOutput with gradient
                biasOutputGradient[k] += partialDerivativeErrorRateWrtBiasOutput;
            }

            printf("\n");

            // update weights and biases
            for (int k = 0; k < numOutputs; k++) {
                for (int l = 0; l < numHidden; l++) {
                    printf("\toutputWeightsForEachHiddenNeuronGradient (outputNum: %d, hiddenNum: %d): %f\n", k, l, outputWeightsForEachHiddenNeuronGradient[k][l]);
                    outputWeightsForEachHiddenNeuron[k][l] += outputWeightsForEachHiddenNeuronGradient[k][l];
                }
                printf("\tbiasOutputGradient (outputNum: %d): %f\n", k, biasOutputGradient[k]);
                biasOutput[k] += biasOutputGradient[k];
            }
            printf("\n");

            for (int k = 0; k < numHidden; k++) {
                for (int l = 0; l < numInputs; l++) {
                    printf("\thiddenWeightsForEachInputNeuronGradient (inputNum: %d, hiddenNum: %d): %f\n", l, k, hiddenWeightsForEachInputNeuronGradient[k][l]);
                    hiddenWeightsForEachInputNeuron[k][l] += hiddenWeightsForEachInputNeuronGradient[k][l];
                }
                printf("\tbiasHiddenGradient (hiddenNum: %d): %f\n", k, biasHiddenGradient[k]);
                biasHidden[k] += biasHiddenGradient[k];
            }
        }

        printf("\n");
        printf("epoch: %d errorRate: %f\n", i, totalError);
    }

    // simple forward propagation after training
    printf("\n");
    printf("STEP_1: Calculate accuracy after training:\n");
    printf("---------------\n");
    int correctMatchedAfterTraining = 0;
    for (int i = 0; i < totalTrainingData; i++) {

        // print input data
        printf("training data: %d\n", i);
        for (int j = 0; j < numInputs; j++) {
            printf("\t%d: %f \n", j, trainingInputs[i][j]);
        }
        printf("\n");

        // calculate hidden outputs
        double computedHiddenNeuronOutput[numHidden];

        for (int j = 0; j < numHidden; j++) {
            double sum = 0;
            for (int k = 0; k < numInputs; k++) {
                sum += trainingInputs[i][k] * hiddenWeightsForEachInputNeuron[j][k];
            }
            sum += biasHidden[j];
            computedHiddenNeuronOutput[j] = sigmoid(sum);
        }

        // calculate output outputs
        double computedOutputNeuronOutput[numOutputs];

        for (int j = 0; j < numOutputs; j++) {
            double sum = 0;
            for (int k = 0; k < numHidden; k++) {
                sum += computedHiddenNeuronOutput[k] * outputWeightsForEachHiddenNeuron[j][k];
            }
            sum += biasOutput[j];
            computedOutputNeuronOutput[j] = sigmoid(sum);
        }

        // find max value in computedOutputNeuronOutput
        double max = 0;
        int maxIndex = 0;
        for (int j = 0; j < numOutputs; j++) {

            // printf("computedOutputNeuronOutput[%d] %f", j, computedOutputNeuronOutput[j]);

            if (computedOutputNeuronOutput[j] > max) {
                max = computedOutputNeuronOutput[j];
                maxIndex = j;
            }
        }

        // print computed output
        printf("\tcomputed output: \n");
        for (int j = 0; j < numOutputs; j++) {
            printf("\t\t%d: %f \n", j, computedOutputNeuronOutput[j]);
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

        printf("\n");

        double correctOutput[numOutputs];

        // get correct output from correctTrainingOutputs
        for (int j = 0; j < numOutputs; j++) {
            correctOutput[j] = correctTrainingOutputs[i][j];
        }

        printf("\tcorrect output: \n");
        for (int j = 0; j < numOutputs; j++) {
            printf("\t\t%d: %f \n", j, correctOutput[j]);
        }

        printf("\t------------\n\n");
        
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
            printf("\tCorrect\n");
            correctMatchedAfterTraining++;
        } else {
            printf("\tIncorrect\n");
        }

        printf("\n");
    }

    // print accuracy
    printf("after training accuracy: %d%%\n", (int)((double) correctMatchedAfterTraining * 100 / totalTrainingData));
    printf("\n");
}