#include <iostream>
#include "trainData.h"
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iterator>
#include <algorithm>
using namespace std;
Traindata::Traindata(char fileName[])
{
    string line; // чтение из файла
    char * cstr = new char [line.length()+1];
    ifstream readFile;
    readFile.open(fileName, ios_base::out);
    if (!readFile.is_open())
    std::cout << "ERROR! " << fileName << " isnt open!" << endl;
    else
    while (getline(readFile,line))
    {
        strcpy (cstr, line.c_str());
        char * chek = std::strtok (cstr,"\t");
        chek = strtok(NULL,"\t");
        int tempTime = std::atoi(cstr);
        int tempCheck = std::atoi(chek);
        inputTime.push_back(tempTime);      // запись в вектор времени
        avgChek.push_back(tempCheck);       // зпись в вектор среднего чека
    }
    inputDataArraySize = avgChek.size();
    float* temp = new float[inputDataArraySize];
    for (int i = 0; i < inputDataArraySize; i++)
    {
        temp[i] = avgChek[i]*0.05;
        profit.push_back(temp[i]);          // запись в вектор прибыли
    }
    readFile.close();
}
Traindata::~Traindata()
{
}
void Traindata::normData()
{
    float* tempProfit = new float[inputDataArraySize];
    float* tempAvgchek = new float[inputDataArraySize];
    float* tempTime = new float[inputDataArraySize];
    float* tempTimeCopy = new float[inputDataArraySize];
    for (int i = 1; i < inputDataArraySize; i++)
    {
        tempTime[i] = i;
    }
    std::vector<float>:: const_iterator avgchekMax = max_element(avgChek.begin(), avgChek.end());
    std::vector<float>:: const_iterator avgchekMin = min_element(avgChek.begin(), avgChek.end());
    std::vector<float>:: const_iterator profitMax = max_element(profit.begin(), profit.end());
    std::vector<float>:: const_iterator profitMin = min_element(profit.begin(), profit.end());
    for (int i = 0; i < inputDataArraySize; i++)
    {
        tempProfit[i] = (profit[i] - *profitMin)/(*profitMax - *profitMin);
        tempAvgchek[i] = (avgChek[i] - *avgchekMin)/(*avgchekMax - *avgchekMin);
        tempTimeCopy[i] = (tempTime[i] - tempTime[0])/(tempTime[inputDataArraySize]-tempTime[0]);
        profitNorm.push_back(tempProfit[i]);
        avgchekNorm.push_back(tempAvgchek[i]);
        timeNorm.push_back(tempTimeCopy[i]);
    }
    delete [] tempAvgchek;
    delete [] tempProfit;
    delete [] tempTime;
    delete [] tempTimeCopy;
}
/* Матрица сврбодной переменной*/
void Traindata::freeVar()
{
}
    /* Транспонированиая Матрица сврбодной переменной*/
void Traindata::freeVarT()
{
}
    /* Перемножение матриц (tranA*A) */
void Traindata::multMatrx()
{
}
void Traindata::determ()
{
}
void Traindata::linReg()
{
}
