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
    void linear_Regresion();        // Функция осущ линейную регрессию
    void normData();                // Нормализация данных
    void freeVar();                 // Матрица свободных переменных
    void freeVarT();                // Транспонированая матрица свободных переменных
    void multMatrx();               // Перемножение матриц
    void inverceMatrix();           // Обратная матрицы
/**
    // freeVar_number       - number of " Omega's "
    // rows, cols           - number of ROWS and COLUMNS for GENERAL DATA ARRAY.
    //                        * rows - total number of each FREE VARIABLE
    //                        * cols - total number of FREE VARIABLES
    // input_time           - vector of time values. Size [rows]
    // time_norm            - normalized vector of time values. Size [rows].
    // average_check        - vector of average check. Size [rows].
    // average_check_norm   - normalized vector of average check. Size [rows].
    // customer             - vector with number of customers. Size [rows].
    // customer_norm        - normalized vector with number of customers. Size [rows].
    // mfreeVar             - matrix of FREE VARIABLES. Size [rows][freeVar_number+1].
    // mfreeVarT            - transposed matrix of FREE VARIABLES. Size [freeVar_number+1][rows].
    // matrixMult           - result of |mfreeVarT*mfreeVar|. Size [freeVar_number+1][freeVar_number+1].
    // inverce_Matrix       - inversed matrixMult.
    // dataTable            - GENERAL DATA ARRAY. Size [rows][cols].
    // weights              - vector of weights [Omega1, Omega2,...OmegaN]. Size [freeVar_number+1].
*/
private:
    std::vector <float> input_time;            // Вектор с данными времени
    std::vector <float> time_norm;
    std::vector <float> average_check;            // Вектор с данными о среднем чеке
    std::vector <float> average_check_norm;
    std::vector <float> profit;             // Вектор с данными о чистой прибыли
    std::vector <float> profitNorm;
    std::vector <float> customer;            // вектор покупателей
    std::vector <float> customer_norm;
    float   **mfreeVar,                     // Матрица свободной переменной
            **mfreeVarT,                    // транспонированая Матрица свободной переменной
            **matrixMult,                   // результирующая матрица mfreeVarT*mfreeVar
            **inverce_Matrix,               // обратная матрица matrixMult
             *weights,                      // искомый вектор параметров ОМЕГА
            **dataTable;                    // массив со вводными данными
    int     rows,
            cols,
            freeVar_number;                 // количество свободных переменных
    size_t inputDataArraySize;        // Размер массива обрабатываемых данных. Нужен для работы счетчиков
};
#endif // TRAINDATA_H_INCLUDED
