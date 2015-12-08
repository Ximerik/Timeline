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
    std::cout <<	 "ERROR! " << fileName << " isnt open!" << endl;
    else
	    while (getline(readFile,line))
	    {
	        strcpy (cstr, line.c_str());
	        char * chek = std::strtok (cstr,"\t");
	        chek = strtok(NULL,"\t");
	        size_t tempTime = std::atoi(cstr);
	        size_t tempCheck = std::atoi(chek);
	        inputTime.push_back(tempTime);      // запись в вектор времени
	        avgChek.push_back(tempCheck);       // зпись в вектор среднего чека
	    }
    inputDataArraySize = avgChek.size();
    float* temp = new float[inputDataArraySize];
    for (size_t i = 0; i < inputDataArraySize; i++)
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
    for (size_t i = 0; i < inputDataArraySize; i++)
    {
        tempTime[i] = 0;
        timeNorm.push_back(i);
    }
    std::vector<float>:: const_iterator avgchekMax = max_element(avgChek.begin(), avgChek.end());
    std::vector<float>:: const_iterator avgchekMin = min_element(avgChek.begin(), avgChek.end());
    std::vector<float>:: const_iterator profitMax = max_element(profit.begin(), profit.end());
    std::vector<float>:: const_iterator profitMin = min_element(profit.begin(), profit.end());
    std::vector<float>:: const_iterator timeNormMin = min_element(timeNorm.begin(), timeNorm.end());
    std::vector<float>:: const_iterator timeNormMax = max_element(timeNorm.begin(), timeNorm.end());
    for (size_t i = 0; i < inputDataArraySize; i++)
    {
        tempProfit[i] = (profit[i] - *profitMin)/(*profitMax - *profitMin);
        tempAvgchek[i] = (avgChek[i] - *avgchekMin)/(*avgchekMax - *avgchekMin);
        tempTime[i] = (timeNorm[i] - *timeNormMin)/(*timeNormMax - *timeNormMin);
        profitNorm.push_back(tempProfit[i]);
        avgchekNorm.push_back(tempAvgchek[i]);
        timeNorm[i] = tempTime[i];
    }
    delete [] tempAvgchek;
    delete [] tempProfit;
    delete [] tempTime;
}
/* Матрица сврбодной переменной*/
void Traindata::freeVar()
{
    mfreeVar = new float*[inputDataArraySize];
    for (size_t i = 0; i < inputDataArraySize; i++)
        mfreeVar[i] = new float[2];
    for (size_t i = 0; i < inputDataArraySize; i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
            mfreeVar[i][j] = 1;
        }
    }
    for (size_t i = 0; i < inputDataArraySize; i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
            mfreeVar[i][1] = timeNorm[i];
        }
    }
}
    /* Транспонированиая Матрица сврбодной переменной*/
void Traindata::freeVarT()
{
    mfreeVarT = new float*[2];
    for (size_t i = 0; i < 2; i++)
        mfreeVarT[i] = new float[inputDataArraySize];
    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < inputDataArraySize; j++)
        {
            mfreeVarT[i][j] = 1;
        }
    }
    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < inputDataArraySize; j++)
        {
            mfreeVarT[1][j] = timeNorm[j];
        }
    }

}
    /* Перемножение матриц (tranA*A) */
void Traindata::multMatrx()
{
    float a,b;
    matrixMult = new float*[2];
    for (size_t i = 0; i < 2; i++)
        matrixMult[i] = new float[2];
    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
	        matrixMult[i][j] = 0;
	            for (size_t k = 0; k < 2; k++)
	                {
	                    a = mfreeVarT[i][k];
	                    b = mfreeVar[k][j];
	                }
        	matrixMult[i][j] = a*b;
        	cout << matrixMult[i][j] << " ";
        }
    }
}
void Traindata::determ()
{
    float a,b,c,d;
    if (sizeof(matrixMult)==4)
    {
        for (size_t i = 0; i < 2; i++)
            {
                for (size_t j = 0; j < 2; j++)
                {
                    a = matrixMult[0][0];
                    b = matrixMult[1][1];
                    a = matrixMult[0][1];
                    a = matrixMult[1][0];
                    determinant = a*b-c*d;
                }
            }
    }
    else
    {
        cout << " THIS fucking bullshit! ";
    }
}
void Traindata::linReg()
{
}
