#include <iostream>
#include "trainData.h"
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <sstream>

using namespace std;
Traindata::Traindata(char fileName[])
{
    freeVar_number = 2;
    string line, buf;
    std::vector<string> tokens;
    std::vector<int> tempCols;
        rows = 0,
        cols = 0;

    ifstream readFile;
    readFile.open(fileName, ios_base::out);
        if (!readFile.is_open())
            std::cout << "ERROR! " << fileName << " isnt open!" << endl;
        else
            while(getline(readFile,line))
            {
                stringstream ss(line);
                while (ss >> buf)
                {
                    tokens.push_back(buf);
                }
                rows++;
                tempCols.push_back(std::count(line.begin(),line.end(),'\t')+1);
                cols = *std::max_element(tempCols.begin(), tempCols.end());
            }
    readFile.close();

    int s = tokens.size();
    float *generalData = new float[s];
    for (int i = 0; i < s; i++)
        {
            generalData[i] = 0;
            generalData[i] = atof(tokens[i].c_str());
        }

    dataTable = new float*[rows];
    for (int i = 0; i < rows; i++)
        dataTable[i] = new float[cols];

    int k = 0;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            dataTable[i][j] = 0;
            dataTable[i][j] = generalData[k];
            k++;
        }

    mfreeVar = new float*[rows];
    for (int i = 0; i < rows; i++)
        mfreeVar[i] = new float[freeVar_number+1];
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < freeVar_number+1; j++)
        {
            mfreeVar[i][j] = 0;
        }

    mfreeVarT = new float*[freeVar_number+1];
    for (int i = 0; i < freeVar_number+1; i++)
        mfreeVarT[i] = new float[rows];
    for (int i = 0; i < freeVar_number+1; i++)
        for (int j = 0; j < rows; j++)
        {
            mfreeVar[i][j] = 0;
        }

    matrixMult = new float*[freeVar_number+1];
    for (int i = 0; i < freeVar_number+1; i++)
        matrixMult [i] = new float[freeVar_number+1];
    for (int i = 0; i < freeVar_number+1; i++)
        for (int j = 0; j < freeVar_number+1; j++)
        {
            matrixMult[i][j] = 0;
        }

    delete [] generalData;
    generalData = nullptr;
}
Traindata::~Traindata()
{
}

void Traindata::dataVectors()
{
    float   time_Element,
            customer_Element;
    for (int i = 0;  i < rows; i++)
    {
        customer_Element = dataTable[i][2];
        time_Element = i;
        customer.push_back(customer_Element);
        inputTime.push_back(time_Element);
    }

}

void Traindata::normData()
{
    std::vector<float>:: const_iterator customer_Max = max_element(customer.begin(), customer.end());
    std::vector<float>:: const_iterator customer_Min = min_element(customer.begin(), customer.end());
    std::vector<float>:: const_iterator inputTime_Max = max_element(inputTime.begin(), inputTime.end());
    std::vector<float>:: const_iterator inputTime_Min = min_element(inputTime.begin(), inputTime.end());

    for (int i = 0;  i < rows; i++)
    {
        timeNorm        .push_back((inputTime[i] - *inputTime_Min)/(*inputTime_Max-*inputTime_Min));
        customerNorm    .push_back((customer[i] - *customer_Min)/(*customer_Max -*customer_Min));
    }
}

void Traindata::freeVar()
{
    for (int i = 0; i < rows; i++)
    {
        mfreeVar[i][0] = 1;
        mfreeVar[i][1] = timeNorm[i];
        mfreeVar[i][2] = customerNorm[i];
    }
}

void Traindata::freeVarT()
{
    for (int i = 0; i < rows; i++)
    {
        mfreeVarT[0][i] = 1;
        mfreeVarT[1][i] = timeNorm[i];
        mfreeVarT[2][i] = customerNorm[i];
    }
}

void Traindata::multMatrx()
{
    int k = 0;
    for (int i = 0; i < freeVar_number+1; i++)
    {
        ++k;
        for (int j = 0; j < freeVar_number+1; j++)
        {
            matrixMult[i][j] += mfreeVarT[i][k]*mfreeVar[k][j];
            cout << matrixMult[i][j] << "\t";
        }
        cout << endl;
    }

}
