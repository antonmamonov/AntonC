#include <stdio.h>
#include <math.h>
#include <stdbool.h>

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
    const int totalTrainingData = 2;

    double trainingInputs[totalTrainingData][numInputs] = {
        {0, 0},
        {0, 1},
        // {1, 0},
        // {1, 1}
    };

    double correctTrainingOutputs[totalTrainingData][numOutputs] = {
        {1, 0},
        {0, 1},
        // {0, 1},
        // {1, 0}
    };

    // initialize all values to zero instead of randomizing so the output is the same every time

    // initialize hidden weights and bias to defaultWeight

    double defaultWeight = 0.420;
    for (int i = 0; i < numHidden; i++) {
        biasHidden[i] = defaultWeight;
        for (int j = 0; j < numInputs; j++) {
            weightsInputToHidden[j][i] = defaultWeight;
        }
    }

    // initialize output weights and bias to 0.5
    for (int i = 0; i < numOutputs; i++) {
        biasOutput[i] = defaultWeight;
        for (int j = 0; j < numHidden; j++) {
            weightsHiddenToOutputs[j][i] = defaultWeight;
        }
    }

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

        double weightsHiddenToOutputsDelta[numHidden][numOutputs];

        // initialize weightsHiddenToOutputsDelta to 0
        for (int j = 0; j < numHidden; j++) {
            for (int k = 0; k < numOutputs; k++) {
                weightsHiddenToOutputsDelta[j][k] = 0;
            }
        }

        double biasOutputDelta[numOutputs];

        // initialize biasOutputDelta to 0
        for (int j = 0; j < numOutputs; j++) {
            biasOutputDelta[j] = 0;
        }

        double totalError = 0;
        for (int j = 0; j < totalTrainingData; j++) {

            printf("\n");
            printf("training data: %d\n", j);
            printf("\n");

            // calculate hidden layer
            double computedHiddenLayer[numHidden];

            for (int k = 0; k < numHidden; k++) {
                double sum = 0;
                for (int l = 0; l < numInputs; l++) {
                    sum += trainingInputs[j][l] * weightsInputToHidden[l][k];

                    // printf("\ttrainingInputs[%d][%d]: %f", j, l, trainingInputs[j][l]);
                    // printf(" * weightsInputToHidden[%d][%d]: %f", l, k, weightsInputToHidden[l][k]);
                    // printf(" = %f", trainingInputs[j][l] * weightsInputToHidden[l][k]);
                    // printf("\n");
                }
                // printf("\tbiasHidden: %f\n", biasHidden[k]);
                // printf("\tsum: %f\n", sum);
                sum += biasHidden[k];
                // printf("\tsumAfterBias: %f\n", sum);
                computedHiddenLayer[k] = sigmoid(sum);
                // printf("\tcomputedHiddenLayer with sigmoid[%d]: %f \n", k, computedHiddenLayer[k]);
            }

            printf("\n");

            // calculate output layer
            double computedOutputLayer[numOutputs];

            for (int k = 0; k < numOutputs; k++) {
                double sum = 0;
                for (int l = 0; l < numHidden; l++) {
                    sum += computedHiddenLayer[l] * weightsHiddenToOutputs[l][k];

                    printf("\tcomputedHiddenLayer[%d]: %f", l, computedHiddenLayer[l]);
                    printf(" * weightsHiddenToOutputs[%d][%d]: %f", l, k, weightsHiddenToOutputs[l][k]);
                    printf(" = %f", computedHiddenLayer[l] * weightsHiddenToOutputs[l][k]);
                    printf("\n");
                }
                printf("\tsum: %f\n", sum);
                printf("\tbiasOutput: %f\n", biasOutput[k]);
                printf("\tbiasOutputPointer: %p\n", &biasOutput[k]);
                sum += biasOutput[k];
                printf("\tOutputLayerSumAfterBias: %f\n", sum);
                computedOutputLayer[k] = sigmoid(sum);
                printf("\tsigmoid: %f\n", sigmoid(sum));
                printf("\tcomputedOutputLayer with sigmoid[%d]: %f\n", k, computedHiddenLayer[k]);
                printf("\n");
            }
            
            // print computed output
            printf("\tcomputed output for training data %d: \n", j);
            for (int k = 0; k < numOutputs; k++) {
                printf("\t\t%d: %f \n", k, computedOutputLayer[k]);
            }

            // ok, now let's do backpropagation!

            // calculate output error
            double error[numOutputs];

            for (int k = 0; k < numOutputs; k++) {

                double errorRate = (correctTrainingOutputs[j][k] - computedOutputLayer[k]);

                printf("\n");
                printf("\t\tcorrectTrainingOutputs (outputNum: %d): %f\n", k, correctTrainingOutputs[j][k]);
                printf("\t\tcomputedOutputLayer (outputNum: %d): %f\n", k, computedOutputLayer[k]);
                printf("\t\terrorRate %f\n\n", errorRate);

                error[k] = errorRate;

                totalError += fabs(errorRate);
            }

            // calculate output delta
            double deltaOutput[numOutputs];

            for (int k = 0; k < numOutputs; k++) {
                double sigmoidDerivativeOutput = sigmoidDerivative(computedOutputLayer[k]);
                deltaOutput[k] = error[k] * sigmoidDerivativeOutput;

                printf("\t\tsigmoidDerivativeOutput (outputNum: %d): %f\n", k, sigmoidDerivativeOutput);
                printf("\t\tdeltaOutput (outputNum: %d): %f\n", k, deltaOutput[k]);
            }

            printf("\n");

            // calculate delta bias output
            double deltaBiasOutput[numOutputs];

            for (int k = 0; k < numOutputs; k++) {
                double sigmoidDerivativeOutputBias = sigmoidDerivative(computedOutputLayer[k]);
                deltaBiasOutput[k] = error[k] * sigmoidDerivativeOutputBias;

                printf("\t\tdeltaBiasOutput (outputNum: %d): %f\n", k, deltaBiasOutput[k]);
            }

            // update delta weights for hidden to output
            for (int k = 0; k < numHidden; k++) {
                for (int l = 0; l < numOutputs; l++) {
                    weightsHiddenToOutputsDelta[k][l] += learningRate * deltaOutput[l];
                }
            }

            // update delta bias output
            for (int k = 0; k < numOutputs; k++) {
                biasOutputDelta[k] += learningRate * deltaBiasOutput[k];
            }

            // calculate hidden layer error
            // double hiddenLayerError[numHidden];

            // for (int k = 0; k < numHidden; k++) {
            //     double sum = 0;
            //     for (int l = 0; l < numOutputs; l++) {
            //         sum += deltaOutput[l] * weightsHiddenToOutputs[k][l];
            //     }
            //     hiddenLayerError[k] = sum;
            // }

            // // calculate delta hidden
            // double deltaHidden[numHidden];

            // for (int k = 0; k < numHidden; k++) {
            //     deltaHidden[k] = hiddenLayerError[k] * sigmoidDerivative(computedHiddenLayer[k]);
            // }

            // // update weights for input to hidden
            // for (int k = 0; k < numInputs; k++) {
            //     for (int l = 0; l < numHidden; l++) {
            //         // weightsInputToHidden[k][l] += learningRate * deltaHidden[l];
            //     }
            // }
        }

        printf("\n");
        // update weights for hidden to output from delta weights
        for (int j = 0; j < numHidden; j++) {
            for (int k = 0; k < numOutputs; k++) {
                weightsHiddenToOutputs[j][k] += weightsHiddenToOutputsDelta[j][k];

                printf("\tweightsHiddenToOutputsDelta (hiddenNum: %d, outputNum: %d): %f\n", j, k, weightsHiddenToOutputsDelta[j][k]);
                printf("\tnew weightsHiddenToOutputs (hiddenNum: %d, outputNum: %d): %f\n", j, k, weightsHiddenToOutputs[j][k]);
            }
        }

        printf("\n");
        // update bias output from delta bias output
        for (int j = 0; j < numOutputs; j++) {
            biasOutput[j] += biasOutputDelta[j];

            printf("\tbiasOutputDelta (outputNum: %d): %f\n", j, biasOutputDelta[j]);
            printf("\tnew biasOutput (outputNum: %d): %f\n", j, biasOutput[j]);
        }

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