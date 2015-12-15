#ifndef TRAINDATA_H_INCLUDED
#define TRAINDATA_H_INCLUDED
#include <vector>

using namespace std;

class Traindata
{


public:
    Traindata(char*);
    ~Traindata();
    void dataVectors();             // распределение данных по векторам
    void linReg();                  // Функция осущ линейную регрессию
    void normData();                // Нормализация данных
    void freeVar();                 // Матрица свободных переменных
    void freeVarT();                // Транспонированая матрица свободных переменных
    void multMatrx();               // Перемножение матриц
    void invrMatrx();               // Обратная матрицы

private:
    std::vector <float> inputTime;            // Вектор с данными времени
    std::vector <float> timeNorm;
    std::vector <float> avgChek;            // Вектор с данными о среднем чеке
    std::vector <float> avgchekNorm;
    std::vector <float> profit;             // Вектор с данными о чистой прибыли
    std::vector <float> profitNorm;
    std::vector <float> customer;            // вектор покупателей
    std::vector <float> customerNorm;
    float   **mfreeVar,                     // Матрица свободной переменной
            **mfreeVarT,                    // транспонированая Матрица свободной переменной
            **matrixMult,                   // результирующая матрица mfreeVarT*mfreeVar
            **invMatrix,                    // обратная матрица matrixMult
            **dataTable;                    // массив со вводными данными
    int     rows,
            cols,
            freeVar_number;                 // количество свободных переменных
    size_t inputDataArraySize;        // Размер массива обрабатываемых данных. Нужен для работы счетчиков
};


#endif // TRAINDATA_H_INCLUDED
