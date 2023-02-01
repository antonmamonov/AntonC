# Neural Network in C

## Basic XOR (ie, "Hello World" of neural networks)

```console
$ gcc basicXOR.c && ./a.out
```

## Training

<img src='images/errorRateVsEpoch.png'/>

Error Rate - Vertical/Y Axis

Epoch step - Horizontal/X Axis

## Theoretical Intuition

Error Rate Function

<img src='images/errorRateFunction.png'>

Sigmoid Activation function

<img src='images/sigmoid.png'>

Partial Derivative of error rate with respect to x

<img src='images/partialDwrtX.png'>

Sigmoid activation with weight (x) given sibling neuron outputs + bias (b) and direct input (c)

<img src='images/sigmoidWithBiasAndInputWeight.png'>

Partial Derivative of sigmoid activation with respect to weight (x) given sibling neuron outputs + bias (b) and direct input (c)

<img src='images/partialDwrtXsigmoidWithBiasAndInputWeight.png'>

Chain Rule

<img src='images/chainRule.png'>
