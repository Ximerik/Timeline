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
    //obj.invrMatrx();
    return 0;
}
