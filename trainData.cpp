#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <sstream>
#include "trainData.h"
#include "linalg.h"
#include "ap.h"
#include <valarray>

using namespace std;

Traindata::Traindata(char fileName[])
{
    /** ==The definition of variables==
    //
    // Global variables:
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

    freeVar_number = 2;
    string line, buf; //
    std::vector<string> tokens;
    std::vector<int> tempCols;
        rows = 0,
        cols = 0;
    // Write data for source file to GENERAL DATA ARRAY
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

    // The definition of global variables
    mfreeVar = new float*[rows];
    for (int i = 0; i < rows; i++)
        mfreeVar[i] = new float[freeVar_number+1];
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < freeVar_number+1; j++)
            mfreeVar[i][j] = 0;

    mfreeVarT = new float*[freeVar_number+1];
    for (int i = 0; i < freeVar_number+1; i++)
        mfreeVarT[i] = new float[rows];
    for (int i = 0; i < freeVar_number+1; i++)
        for (int j = 0; j < rows; j++)
            mfreeVar[i][j] = 0;

    matrixMult = new float*[freeVar_number+1];
    for (int i = 0; i < freeVar_number+1; i++)
        matrixMult [i] = new float[freeVar_number+1];
    for (int i = 0; i < freeVar_number+1; i++)
        for (int j = 0; j < freeVar_number+1; j++)
            matrixMult[i][j] = 0;

    inverce_Matrix = new float*[freeVar_number+1];
    for(int i = 0; i < freeVar_number+1; i++)
        inverce_Matrix[i] = new float [freeVar_number+1];
    for(int i = 0; i < freeVar_number+1; i++)
        for(int j = 0; j < freeVar_number+1; j++)
            inverce_Matrix[i][j] = 0;

    weights = new float [freeVar_number+1];
    for(int i = 0; i < freeVar_number+1; i++)
        weights[i] = 0;
    // Free memory
    delete [] generalData;
    generalData = nullptr;
}
Traindata::~Traindata()
{
}

void Traindata::dataVectors()
{
    /** Write FREE VARIABLES from GENERAL DATA ARRAY to independent vectors */
    float   time_Element,
            customer_Element,
            average_Check;
    for (int i = 0;  i < rows; i++)
    {
        customer_Element = dataTable[i][2];
        average_Check = dataTable[i][1];
        time_Element = dataTable[i][0];
        customer.push_back(customer_Element);
        input_time.push_back(time_Element);
        average_check.push_back(average_Check);
    }

}

void Traindata::normData()
{
    /** Define normalized values of FREE VARIABLES vectors

        VNi = (Vi - Vmin)/(Vmax - Vmin)

        VNi     - normalized value;
        Vi      - current element vector;
        Vmin    - minimum element of vector;
        Vmax    - maximum element of vector.
    */
    std::vector<float>:: const_iterator customer_Max = max_element(customer.begin(), customer.end());
    std::vector<float>:: const_iterator customer_Min = min_element(customer.begin(), customer.end());
    std::vector<float>:: const_iterator inputTime_Max = max_element(input_time.begin(), input_time.end());
    std::vector<float>:: const_iterator inputTime_Min = min_element(input_time.begin(), input_time.end());

    for (int i = 0;  i < rows; i++)
    {
        time_norm        .push_back((input_time[i] - *inputTime_Min)/(*inputTime_Max-*inputTime_Min));
        customer_norm    .push_back((customer[i] - *customer_Min)/(*customer_Max -*customer_Min));
    }
}

void Traindata::freeVar()
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < freeVar_number+1; j++)
        {
            mfreeVar[i][0] = 1;
            mfreeVar[i][1] = input_time[i];
            mfreeVar[i][2] = customer[i];
            //cout << mfreeVar[i][j] << " ";
        }   //cout << endl;
    }

}

void Traindata::freeVarT()
/*************************************************************************
Cache-oblivous real "copy-and-transpose"

Input parameters:
    M   -   number of rows
    N   -   number of columns
    A   -   source matrix, MxN submatrix is copied and transposed
    IA  -   submatrix offset (row index)
    JA  -   submatrix offset (column index)
    B   -   destination matrix, must be large enough to store result
    IB  -   submatrix offset (row index)
    JB  -   submatrix offset (column index)
*************************************************************************/
{
    alglib::ae_int_t m=rows,
                     n=freeVar_number+1,
                     ia=0,
                     ja=0,
                     ib=0,
                     jb=0;
    alglib::real_2d_array matrix_one,
                          matrix_two;
    matrix_one.setlength(rows,freeVar_number+1);
    matrix_two.setlength(freeVar_number+1,rows);
        for (int i = 0; i < rows; i++)
                for (int j = 0; j < freeVar_number+1; j++)
                    matrix_one(i,j)=mfreeVar[i][j];
    alglib::rmatrixtranspose(m,n,matrix_one,
                             ia,ja,matrix_two,
                             ib,jb);
    for (int j = 0; j < freeVar_number+1; j++)
    {
        for (int i = 0; i < rows; i++)
        {
            mfreeVarT[j][i]=matrix_two(j,i);
        }

    }



}

void Traindata::multMatrx()
    /** FUNCTION MANUAL
    // rmatrixgemm() function allows us to calculate matrix product C:=A*B or
    // to perform more general operation, C:=alpha*op1(A)*op2(B)+beta*C,
    // where A, B, C are rectangular matrices, op(X) can be X or X^T,
    // alpha and beta are scalars.
    //
    // This function:
    // * can apply transposition and/or multiplication by scalar to operands
    // * can use arbitrary part of matrices A/B (given by submatrix offset)
    // * can store result into arbitrary part of C
    // * for performance reasons requires C to be preallocated
    //
    // Parameters of this function are:
    // * M, N, K            -   sizes of op1(A) (which is MxK), op2(B) (which
    //                          is KxN) and C (which is MxN)
    // * Alpha              -   coefficient before A*B
    // * A, IA, JA          -   matrix A and offset of the submatrix
    // * OpTypeA            -   transformation type:
    //                          0 - no transformation
    //                          1 - transposition
    // * B, IB, JB          -   matrix B and offset of the submatrix
    // * OpTypeB            -   transformation type:
    //                          0 - no transformation
    //                          1 - transposition
    // * Beta               -   coefficient before C
    // * C, IC, JC          -   preallocated matrix C and offset of the submatrix
    //
    // Below we perform simple product C:=A*B (alpha=1, beta=0)
    //
    // IMPORTANT: this function works with preallocated C, which must be large
    //            enough to store multiplication result.
    //
    // -- ALGLIB --
    // Copyright 2005-2010 by Bochkanov Sergey
    // (c) ALGLIB http://alglib.sources.ru
    //
    // This library is used under license GPL 2+,
    // GPL license version 2 or at your option any later version.
    // A copy of the GNU General Public License is available
    // at http://www.fsf.org/licensing/licenses
    //
    */
{
    // set up temporary arrays for calculations
    alglib::real_2d_array matrix_one;
    alglib::real_2d_array matrix_two;
    alglib::real_2d_array result_matrix;
    matrix_one.setlength(freeVar_number+1,rows);
    matrix_two.setlength(rows,freeVar_number+1);
    result_matrix.setlength(freeVar_number+1,freeVar_number+1);
    for (int i = 0; i < freeVar_number+1;i++)
        for (int j = 0; j < rows; j++)
        {
            matrix_one(i,j) = mfreeVarT[i][j];
            matrix_two(j,i) = mfreeVar[j][i];
        }
    // Input parameters for calculation function
    alglib::ae_int_t m = freeVar_number+1,
                     n = freeVar_number+1,
                     k = rows,
                     ia = 0,
                     ib = 0,
                     ic = 0,
                     ja = 0,
                     jb = 0,
                     jc = 0,
                     optypea = 0,
                     optypeb = 0;
    double alpha = 1.0,
           beta = 0;
    // Calculation function
    alglib::rmatrixgemm(m,n,k,
                        alpha, matrix_one,    ia,ja,optypea,
                               matrix_two,    ib,jb,optypeb,
                         beta, result_matrix, ic,jc);
    for (int i = 0; i < freeVar_number+1;i++)
        for (int j = 0; j < freeVar_number+1; j++)
            matrixMult[i][j]=result_matrix(i,j);
}
void Traindata::inverceMatrix()
    /** FUNCTION MANUAL
    // Input parameters:
    // A       -   matrix.
    // N       -   size of matrix A (optional) :
    //             * if given, only principal NxN submatrix is processed  and
    //               overwritten. other elements are unchanged.
    //             * if not given,  size  is  automatically  determined  from
    //               matrix size (A must be square matrix)
    //
    // Output parameters:
    // Info    -   return code, same as in RMatrixLUInverse
    // Rep     -   solver report, same as in RMatrixLUInverse
    // A       -   inverse of matrix A, same as in RMatrixLUInverse
    //
    // Result:
    // True, if the matrix is not singular.
    // False, if the matrix is singular.
    //
    // -- ALGLIB --
    // Copyright 2005-2010 by Bochkanov Sergey
    // (c) ALGLIB http://alglib.sources.ru
    //
    // This library is used under license GPL 2+,
    // GPL license version 2 or at your option any later version.
    // A copy of the GNU General Public License is available
    // at http://www.fsf.org/licensing/licenses
    //
    */
{
    alglib::ae_int_t info;
    alglib::matinvreport rep;
    alglib::real_2d_array inv_matrix; // аватар обратной матрицы
    inv_matrix.setlength(freeVar_number+1,freeVar_number+1);
    for (int i = 0; i < freeVar_number+1; i++)
        for (int j = 0; j < freeVar_number+1; j++)
            inv_matrix(i,j) = matrixMult[i][j];
    alglib::rmatrixinverse(inv_matrix,info,rep);
    for (int i = 0; i < freeVar_number+1; i++)
        for (int j = 0; j < freeVar_number+1; j++)
            inverce_Matrix[i][j] = inv_matrix(i,j);
}
void Traindata::linear_Regresion()
    /** FUNCTION MANUAL
    // rmatrixgemm() function allows us to calculate matrix product C:=A*B or
    // to perform more general operation, C:=alpha*op1(A)*op2(B)+beta*C,
    // where A, B, C are rectangular matrices, op(X) can be X or X^T,
    // alpha and beta are scalars.
    //
    // This function:
    // * can apply transposition and/or multiplication by scalar to operands
    // * can use arbitrary part of matrices A/B (given by submatrix offset)
    // * can store result into arbitrary part of C
    // * for performance reasons requires C to be preallocated
    //
    // Parameters of this function are:
    // * M, N, K            -   sizes of op1(A) (which is MxK), op2(B) (which
    //                          is KxN) and C (which is MxN)
    // * Alpha              -   coefficient before A*B
    // * A, IA, JA          -   matrix A and offset of the submatrix
    // * OpTypeA            -   transformation type:
    //                          0 - no transformation
    //                          1 - transposition
    // * B, IB, JB          -   matrix B and offset of the submatrix
    // * OpTypeB            -   transformation type:
    //                          0 - no transformation
    //                          1 - transposition
    // * Beta               -   coefficient before C
    // * C, IC, JC          -   preallocated matrix C and offset of the submatrix
    //
    // Below we perform simple product C:=A*B (alpha=1, beta=0)
    //
    // IMPORTANT: this function works with preallocated C, which must be large
    //            enough to store multiplication result.
    //
    // -- ALGLIB --
    // Copyright 2005-2010 by Bochkanov Sergey
    // (c) ALGLIB http://alglib.sources.ru
    //
    // This library is used under license GPL 2+,
    // GPL license version 2 or at your option any later version.
    // A copy of the GNU General Public License is available
    // at http://www.fsf.org/licensing/licenses
    //
    */
{
    // Set up temporary arrays for calculations
    float ** resultMatrix = new float * [freeVar_number+1];
    float *  weightMatrix = new float [freeVar_number+1];
    for (int i = 0; i < freeVar_number+1; i++)
    {
        resultMatrix[i] = new float[rows];
        weightMatrix[i] = 0;
        for (int j = 0; j < rows; j++)
            resultMatrix[i][j] = 0;
    }
    alglib::real_2d_array matrix_one;   // аватар inverce_Matrix
    alglib::real_2d_array matrix_two;   // аватар mfreeVarT
    alglib::real_2d_array matrix_three; // аватар average_check
    alglib::real_2d_array weight_matrix;// вектор весовых коефициентов
    alglib::real_2d_array result_matrix;// матрица для промежуточных резудльтатов

    matrix_one.setlength(freeVar_number+1,freeVar_number+1);
    matrix_two.setlength(freeVar_number+1,rows);
    matrix_three.setlength(rows,1);
    result_matrix.setlength(freeVar_number+1,rows);
    weight_matrix.setlength(freeVar_number+1,1);

    for (int i = 0; i < freeVar_number+1; i++)
        for (int j = 0; j < freeVar_number+1; j++)
            matrix_one(i,j) = inverce_Matrix[i][j];
    for (int i = 0; i < freeVar_number+1; i++)
        for (int j = 0; j < rows; j++)
        {
            matrix_two(i,j) = mfreeVarT[i][j];
            result_matrix(i,j) = resultMatrix[i][j];
        }
        for (int i = 0; i < rows; i++)
            matrix_three(i,0) = average_check[i];
        for (int i = 0; i < freeVar_number+1; i++)
            weight_matrix(i,0) = weightMatrix[i];
    //Input parameters for calculation
    //[weights]=[inverce_Matrix]*[mfreeVarT]*[result_matrix]
    alglib::ae_int_t m = freeVar_number+1,
                     k = freeVar_number+1,
                     n = rows,
                     ia = 0,
                     ib = 0,
                     ic = 0,
                     ja = 0,
                     jb = 0,
                     jc = 0,
                     optypea = 0,
                     optypab = 0;
    double alpha = 1.0,
           beta = 0.0;

    alglib::rmatrixgemm(m,n,k,
                        alpha, matrix_one,    ia,ja,optypea,
                               matrix_two,    ib,jb,optypab,
                         beta, result_matrix, ic,jc);
    alglib::rmatrixgemm(m,1,n,
                        alpha, result_matrix, ia,ja,optypea,
                               matrix_three,  ib,jb,optypab,
                         beta, weight_matrix, ic,jc);
    for (int i = 0; i < freeVar_number+1; i++)
        weights[i] = weight_matrix(i,0);
    for (int i = 0; i < freeVar_number+1; i++)
        cout << weights[i] << " ";

    delete [] resultMatrix;
    delete [] weightMatrix;
    resultMatrix = nullptr;
    weightMatrix = nullptr;
}
