#pragma once
#ifndef BASE_H
#define BASE_H
namespace MatrixWork
{
    class Util;
    class Value;
    class Matrix;
    class Solution;
    Matrix operator*(const Matrix &v1, const Matrix &v2);
    Matrix operator+(const Matrix &v1, const Matrix &v2);
    Matrix operator-(const Matrix &v1, const Matrix &v2);
    bool operator==(const Matrix &v1, const Matrix &v2);
    bool operator!=(const Matrix &v1, const Matrix &v2);
    Value operator+(const Value &v1, const Value &v2);
    Value operator-(const Value &v1, const Value &v2);
    Value operator*(const Value &v1, const Value &v2);
    Value operator/(const Value &v1, const Value &v2);
    bool operator==(const Value &v1, const Value &v2);
    bool operator!=(const Value &v1, const Value &v2);
    bool operator>(const Value &v1, const Value &v2);
    bool operator<(const Value &v1, const Value &v2);
}

#endif