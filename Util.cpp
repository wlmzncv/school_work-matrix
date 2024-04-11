#include "Util.h"
#include "Value.h"
#include "Base.h"
#include "Matrix.h"
#include <iostream>
#include <stdexcept>

using std::cout;

namespace MatrixWork{
    const double Util::error = 0.0000000001;
    void Util::testLU(){
        Matrix m(3, 3);
        m[0][0] = 10;
        m[0][1] = -7;
        m[0][2] = 0;
        m[1][0] = -3;
        m[1][1] = 2;
        m[1][2] = 6;
        m[2][0] = 5;
        m[2][1] = -1;
        m[2][2] = 5;
        Matrix L, P, U;
        auto r = m.LU(P, L, U);
        cout << "A=\n"
             << m.toString() << "\n";
        if (r == 0)
        {
            cout << "P=\n"
                 << P.toString() << "\n";
            cout << "L=\n"
                 << L.toString() << "\n";
            cout << "U=\n"
                 << U.toString() << "\n";
            cout << "P*A=\n"
                 << (P * m).toString() << "\n";
            cout << "L*U=\n"
                 << (L * U).toString() << "\n";
        }
        else
        {
            cout << "r=" << r << "\n";
        }
    }
    void Util::testQR(){
        Matrix m(3, 3);
        m[0][0] = 1;
        m[0][1] = 1;
        m[0][2] = 1;
        m[1][0] = 0.001;
        m[1][1] = 0.001;
        m[1][2] = 0;
        m[2][0] = 0.001;
        m[2][1] = 0;
        m[2][2] = 0.001;
        Matrix Q, R;
        auto r = m.QR(Q, R);
        cout << "A=\n"
             << m.toString() << "\n";
        if (r == 0)
        {
            cout << "Q=\n"
                 << Q.toString() << "\n";
            cout << "R=\n"
                 << R.toString() << "\n";
            cout << "Q*R=\n"
                 << (Q * R).toString() << "\n";
        }
        else
        {
            cout << "r=" << r << "\n";
        }
    }
    void Util::testHouseholder(){
        Matrix m(4, 3);
        m[0][0] = -4;
        m[0][1] = 2;
        m[0][2] = -4;
        m[1][0] = -2;
        m[1][1] = -2;
        m[1][2] = 1;
        m[2][0] = 4;
        m[2][1] = -2;
        m[2][2] = 4;
        m[3][0] = 2;
        m[3][1] = -1;
        m[3][2] = 2;
        Matrix P, T;
        auto r = m.Householder(P,T);
        cout << "A=\n"
             << m.toString() << "\n";
        if (r == 0)
        {
            cout << "P=\n"
                 << P.toString() << "\n";
            cout << "T=\n"
                 << T.toString() << "\n";
            cout << "P*A=\n"
                 << (P * m).toString() << "\n";
        }
    }
    void Util::testGivens()
    {
        Matrix m(3, 3);
        m[0][0] = -4;
        m[0][1] = -2;
        m[0][2] = 1;
        m[1][0] = -2;
        m[1][1] = -2;
        m[1][2] = -2;
        m[2][0] = 1;
        m[2][1] = -2;
        m[2][2] = 4;
        Matrix P, T;
        auto r = m.Givens(P, T);
        cout << "A=\n"
             << m.toString() << "\n";
        if (r == 0)
        {
            cout << "P=\n"
                 << P.toString() << "\n";
            cout << "T=\n"
                 << T.toString() << "\n";
            cout << "P*A=\n"
                 << (P * m).toString() << "\n";
        }
    }
    void Util::testURV()
    {
        Matrix m(3, 3);
        m[0][0] = -4;
        m[0][1] = -2;
        m[0][2] = 1;
        m[1][0] = -2;
        m[1][1] = -2;
        m[1][2] = -2;
        m[2][0] = 1;
        m[2][1] = -2;
        m[2][2] = 4;
        Matrix U, R, V;
        int r = 0;
        try {
            r = m.URV(U, R, V);
        }catch (std::exception& e){
            cout << e.what() << "\n";
            exit(-1);
        }
        
        cout << "A=\n"
             << m.toString() << "\n";
        if (r == 0)
        {
            cout << "U=\n"
                 << U.toString() << "\n";
            cout << "R=\n"
                 << R.toString() << "\n";
            cout << "V=\n"
                 << V.toString() << "\n";
            cout << "U*R*(V.T)=\n"
                 << (U*R*V.T()).toString() << "\n";
        }
    }
    bool Util::isAllZero(const vector<Value> &v)
    {
        for (Value i : v)
        {
            if (i.abs().toDouble() > Util::error)
            {
                return false;
            }
        }
        return true;
    }
    Value Util::innerProduct(const vector<Value> &v1, const vector<Value> &v2){
        Value t(0,1);
        if(v1.size()!=v2.size()){
            throw std::runtime_error("innerProduct : vectors are not Compatilbe");
        }else {
            auto n = v1.size();
            for (int i = 0; i < n; ++i){
                t += v1[i] * v2[i];
            }
            return t;
        }
    }
    void Util::chengVector(vector<Value> &v, const Value &a){
        for (int i = 0; i < v.size(); ++i){
            v[i] *= a;
        }
    }
    vector<Value> Util::chengVector(unsigned int size, unsigned int n, const Value &a){
        vector<Value> t(size);
        t[n] = a;
        return t;
    }
    int Util::gcd(int a, int b)
    {
        if (a < 0)
            a = -a;
        if (b < 0)
            b = -b;
        if (a < b)
        {
            swap(a, b);
        }
        if (b == 0)
            return a;
        int r = a % b;
        while (r)
        {
            a = b;
            b = r;
            r = a % b;
        }
        return b;
    }
    unsigned int Util::absmaxInVector(const vector<Value> &v, int from){
        if(v.empty())
            return -1;
        else {
            Value max = v[from].abs();
            int l = from;
            for (long i = from+1; i < v.size(); ++i){
                if(v[i].abs() > max){
                    max = v[i].abs();
                    l = i;
                }
            }
            return l;
        }
    }
    Value Util::getValue(const string &s)
    {
        std::istringstream is(s);
        double a, b;
        char c;
        if (s.find("/") != string::npos)
        {
            is >> a >> c >> b;
            return Value(a, b);
        }
        else
        {
            is >> a;
            return Value(a, 1);
        }
    }
}