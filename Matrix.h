#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include "Value.h"
#include "Base.h"

using std::swap;
using std::vector;

class MatrixWork::Matrix
{
    friend Matrix operator*(const Matrix &v1, const Matrix &v2);
    friend bool operator==(const Matrix &v1, const Matrix &v2);
public:
    // 得到单位矩阵
    static Matrix getI(unsigned int size){
        Matrix t(size, size);
        for (long i = 0; i < size; ++i)
        {
            t[i][i] = Value(1, 1);
        }
        return t;
    }
    // 初等变换矩阵 交换两行(列) 左乘该矩阵是做行变换
    static Matrix getE1(unsigned int size, unsigned int a, unsigned int b);
    // 初等变换矩阵 将某一行(列)乘以一个实数
    static Matrix getE2(unsigned int size, unsigned int a, double d);
    // 初等变换矩阵 将a行乘以一个实数加到b行上
    static Matrix getE3(unsigned int size, unsigned int a, unsigned int b, double d);

private:
    unsigned int x; //行数
    unsigned int y; //列数
    vector<vector<Value>> A;

public:
    Matrix():x(0),y(0){};
    Matrix(unsigned int a, unsigned int b);
    Matrix(const Matrix &v):x(v.x),y(v.y){
        for (auto i = 0; i < x; ++i){
            A.push_back(vector<Value>(v[i]));
        }
    }
    Matrix(Matrix &&v) noexcept : x(v.x),y(v.y),A(v.A){
        v.A = {};
        v.x = v.y = 0;
    }
    Matrix(const vector<Value>& v, int flag = 0);
    Matrix &operator=(const Matrix &v);
    Matrix &operator=(Matrix &&v) noexcept;
    ~Matrix();

    vector<Value>& operator[] (std::size_t n){
        return A[n];
    }
    const vector<Value>& operator[] (std::size_t n) const{
        return A[n];
    }
    Matrix &operator+=(const Matrix &v);
    Matrix &operator-=(const Matrix &v);

    Matrix cheng(double a) const;
    Matrix cheng(const Value &a) const;
    vector<Value> columnAt(std::size_t n, int from = 0) const{
        vector<Value> t;
        for (int i = from; i < x; ++i){
            t.push_back(A[i][n]);
        }
        return t;
    }
    Matrix T() const {
        Matrix t;
        t.x = y;
        t.y = x;
        for (int i = 0; i < y; ++i){
            t.A.push_back(columnAt(i));
        }
        return t;
    }
    Matrix frontRows(unsigned int size) const {
        Matrix t;
        if(size <= x){
            for (int i = 0; i < size; ++i){
                t.A.push_back(A[i]);
                (t.x)++;
            }
        }
        if(t.x>0)
            t.y = y;
        return t;
    }
    unsigned int X() const { return x; }
    unsigned int Y() const { return y; }
    bool isSquare() const { return x == y; }
    bool isZero() const {
        for (int i = 0; i < x; ++i){
            if(!Util::isAllZero(A[i])){
                return false;
            }
        }
        return true;
    }
    bool empty() const { return x == 0 || y == 0; }
    std::string toString();
    unsigned int rank() const;

    // 对矩阵A进行LU分解 A=LU
    // 返回0分解成功 1不能对A分解 2分解失败
    int LU(Matrix &P, Matrix &L, Matrix &U) const;
    // 对矩阵A进行QR分解 A=QR
    // 返回0分解成功 1不能对A分解 2分解失败
    int QR(Matrix &Q, Matrix &R) const;
    // 对矩阵A进行Householder约简 PA=T
    // 返回0分解成功
    int Householder(Matrix &P, Matrix &T) const;
    // 对矩阵A进行Givens约简 PA=T
    // 返回0分解成功
    int Givens(Matrix &P, Matrix &T) const;
    // 对矩阵A进行URV分解 A=URV的转置
    // 返回0分解成功
    int URV(Matrix &U, Matrix &R, Matrix &V) const;

private:
    // T为行阶梯阵 将其中不为0的行作为一组基 然后转为标准正交基 按列存储
    Matrix fun1(const Matrix &T) const;
    // 该算法不稳定
};

#endif