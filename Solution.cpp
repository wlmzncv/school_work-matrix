#include "Solution.h"
#include <sstream>
#include <iostream>

namespace MatrixWork{
    using std::cout;
    void Solution::freshState() {
        if(A.X()==0 || A.Y()==0){
            state = "";
            return;
        }
        r = A.rank();
        std::stringstream ss;
        ss << "A's rank is " << r << " ;\n";
        ss << "Can do these Factorization to A: ";
        if(A.X()==A.Y() && A.X() == r){
            ss << "LU, ";
        }
        if(A.Y() == r){
            ss << "QR, ";
        }
        ss << "URV;\n";
        ss << "Can do Householder or Givens Reduction to A;\n";
        if(A.X() == A.Y()){
            ss << "A has Determinant\n";
        }else {
            ss << "A doesn't have Determinant\n";
        }
        state = ss.str();
    }
    void Solution::getDeterminant() {
        if(A.empty()||A.X()!=A.Y()){
            cout << "A doesn't have Determinant\n";
        }else {
            if(r!=A.X()){
                cout << "Determinant is 0\n";
            }else {
                cout << "Use Givens Reduction to calculate Determinant\n";
                A.Givens(P_G, T_G);
                double t=1;
                for (int i = 0; i < T_G.X(); ++i){
                    t *= T_G[i][i].toDouble();
                }
                cout << "Determinant is " << t << " \n";
            }
        }
    }
    void Solution::takeLU() {
        if((!A.empty())&&r==A.X()&&r==A.Y()){
            A.LU(P, L, U);
            cout << "P " << P.toString() << "L " << L.toString() << "U " << U.toString();
            X = Matrix(A.Y(), 1);
            Matrix XX = X;
            fun1(L, XX, P * b, 1);
            fun1(U, X, XX, 0);
            cout << "The Solution for Ax=b is ";
            cout << "X " << X.toString();
        }else {
            cout << "Cannot do LU Factorization\n";
        }
    }
    void Solution::takeQR() {
        if(!A.empty()&&r==A.Y()){
            A.QR(Q, R);
            cout << "Q " << Q.toString() << "R " << R.toString();
            X = Matrix(A.Y(), 1);
            fun1(R, X, Q.T() * b, 0);
            cout << "The Solution for Ax=b is ";
            cout << "X " << X.toString();
        }else {
            cout << "Cannot do QR Factorization\n";
        }
    }
    void Solution::takeHouseholder(){
        if(!A.empty()){
            A.Householder(P_H, T_H);
            cout << "P " << P_H.toString() << "T " << T_H.toString();
            X = Matrix(A.Y(), 1);
            fun1(T_H, X, P_H * b);
            cout << "The Solution for Ax=b is ";
            cout << "X " << X.toString();
        }else {
            cout << "Matrix A is empty\n";
        }
    }
    void Solution::takeGivens(){
        if (!A.empty())
        {
            A.Givens(P_G, T_G);
            cout << "P " << P_G.toString() << "T " << T_G.toString();
            X = Matrix(A.Y(), 1);
            fun1(T_G, X, P_G * b);
            cout << "The Solution for Ax=b is ";
            cout << "X " << X.toString();
        }
        else
        {
            cout << "Matrix A is empty\n";
        }
    }
    void Solution::takeURV(){
        if (!A.empty())
        {
            try{
                A.URV(W, S, V);
            }catch (std::exception &e){
                cout << e.what() << "\n Sorry! Failed.";
                return;
            }
            cout << "U " << W.toString() << "R " << S.toString() << "V " << V.toString();
        }
        else
        {
            cout << "Matrix A is empty\n";
        }
    }
    // 求解线性方程组 A需为上三角或下三角的方阵 flag=0为上三角矩阵
    void Solution::fun1(const Matrix &A, Matrix &X, const Matrix &B, int flag){
        int l = X.X();
        if(A.isSquare()&&A.X()==l&&B.X()==l){
            if(flag == 0){
                for (int i = l - 1; i >= 0; --i)
                {
                    X[i][0] = B[i][0];
                    for (int j = l - 1; j > i; --j)
                    {
                        X[i][0] -= X[j][0] * A[i][j];
                    }
                    X[i][0] /= A[i][i];
                }
            }else {
                for (int i = 0; i < l; ++i)
                {
                    X[i][0] = B[i][0];
                    for (int j = 0; j < i; ++j)
                    {
                        X[i][0] -= X[j][0] * A[i][j];
                    }
                    X[i][0] /= A[i][i];
                }
            }
        }
    }
    void Solution::fun1(const Matrix &R, Matrix &X, const Matrix &C){
        int l = X.X();
        if(l<=R.X() && l<=C.X()){
            fun1(R.frontRows(l), X, C.frontRows(l), 0);
        }
    }
}