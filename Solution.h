#pragma once
#ifndef SOLUTION_H
#define SOLUTION_H

#include "Base.h"
#include "Value.h"
#include "Matrix.h"

class MatrixWork::Solution {
private:
    string state;
    int r;

public:
    // 系数矩阵A
    Matrix A;
    // 线性方程组的右端 行数与A相同 列数为1
    Matrix b;
    // 线性方程组的解 行数与A的列数相同 列数为1
    Matrix X;
    // LU分解
    Matrix P;
    Matrix L;
    Matrix U;
    // QR分解
    Matrix Q;
    Matrix R;
    // Householder约简得到的P和T
    Matrix P_H;
    Matrix T_H;
    // Givens约简得到的P和T
    Matrix P_G;
    Matrix T_G;
    // URV分解 U用W表示 R用S表示
    Matrix W;
    Matrix S;
    Matrix V;

public:
    Solution(){}
    ~Solution(){}
    const string &getState() const { return state; }
    void freshState();
    void takeLU();
    void takeQR();
    void takeHouseholder();
    void takeGivens();
    void takeURV();
    void getDeterminant();

private:
    // 求解线性方程组 A需为上三角或下三角的方阵 flag=0为上三角矩阵
    void fun1(const Matrix &A, Matrix &X, const Matrix &B, int flag);
    void fun1(const Matrix &R, Matrix &X, const Matrix &C);
};

#endif