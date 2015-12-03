#ifndef TRAINDATA_H_INCLUDED
#define TRAINDATA_H_INCLUDED
#include <vector>

using namespace std;

class Traindata
{


public:
    Traindata(char*);
    ~Traindata();
    void linReg();                  // Функция осущ линейную регрессию
    void normData();                // Нормализация данных
    void freeVar();                 // Матрица свободных переменных
    void freeVarT();                // Транспонированая матрица свободных переменных
    void multMatrx();               // Перемножение матриц
    void determ();                  // Определитель матрицы
    void setData();                 // Функция в которой выполняются операции над данными
    void getData();                 // Функция вывода результата

private:
    vector <int> inputTime;            // Вектор с данными времени
    vector <float> timeNorm;
    vector <float> avgChek;            // Вектор с данными о среднем чеке
    vector <float> avgchekNorm;
    vector <float> profit;             // Вектор с данными о чистой прибыли
    vector <float> profitNorm;
    float *A, **tranA, **matrixMult;
    int     inputDataArraySize;        // Размер массива обрабатываемых данных. Нужен для работы счетчиков
};


#endif // TRAINDATA_H_INCLUDED
