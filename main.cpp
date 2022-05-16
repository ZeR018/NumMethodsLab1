#include <iostream>
#include "TMatrix.h"
#include "MMN.h"
#include <Windows.h>

int main()
{
    setlocale(LC_ALL, "Russian");

    // Параметры 
    int N = 15;
    int M = 15;
    int Max_iter = 10000;           // Максимальное число шагов
    double er = 0.0000005;          // Остановка поточности
    TVector <double> xBorder(2);
    TVector <double> yBorder(2);
    TVector <double> accurancy(2);
    TVector <int> max_iter(2);
    accurancy[0] = er;
    accurancy[1] = er;
    max_iter[0] = Max_iter;
    max_iter[1] = Max_iter;
    xBorder[0] = -1;                // a
    xBorder[1] = 1;                 // b 
    yBorder[0] = -1;                // c
    yBorder[1] = 1;                 // d
    std::cout << "Start solving" << std::endl << std::endl;
    MMN mat(N, M, xBorder, yBorder);
    //mat.VectorNevyazki();
    mat.SolvingTestTask(accurancy, max_iter, 'X');
    system("python graph.py");
}