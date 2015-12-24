/**
// ==Abstract==
//
// Linear regression - method of establishing the relationship between the two variables.
// Below is an example of a program that builds a linear model based on a given sample
// and displays the result on the chart.
//
// Model:
// y = Omega1 + Omega2*Xi + Omega2*Xj+...+OmegaN*Xn + Ei + Ej +...+En
//
// Input parameters:
// Omega1, Omega2,...OmegaN  - weights;
// Xi, Xj...Xn               - free variables: time, number of customers etc.
// Ei...En                   - random value. It is assumed that the random variable is
//                             normally distributed with zero math expectancy and
//                             fixed variance.
// Output
// y - average check.
//
*/

#include <iostream>
#include "trainData.h"

using namespace std;

int main()
{
    char fileName[] = "test.txt";
    Traindata obj(fileName);
    obj.dataVectors();
    obj.normData();
    obj.freeVar();
    obj.freeVarT();
    obj.multMatrx();
    obj.inverceMatrix();
    obj.linear_Regresion();
    return 0;
}
