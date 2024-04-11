#include "Matrix.h"
#include <stdexcept>
#include <iostream>

namespace MatrixWork{
    Matrix::Matrix(unsigned int a, unsigned int b):x(a),y(b){
        for (auto i = 0; i < x; ++i){
            A.push_back(vector<Value>(y));
        }
    }
    Matrix::Matrix(const vector<Value>& v, int flag){
        if(flag==0){ //列向量转为矩阵
            x = v.size();
            y = 1;
            for (auto i = 0; i < x; ++i)
            {
                A.push_back(vector<Value>{v[i]});
            }
        }else{ //行向量转为矩阵
            x = 1;
            y = v.size();
            A.push_back(v);
        }
    }
    Matrix & Matrix::operator=(const Matrix &v)
    {
        if(this == &v){
            return *this;
        }
        x = v.x;
        y = v.y;
        for (auto i : A)
        {
            i.clear();
        }
        A.clear();
        for (auto i = 0; i < x; ++i)
        {
            A.push_back(vector<Value>(v[i]));
        }
        return *this;
    }
    Matrix & Matrix::operator=(Matrix &&v) noexcept{
        if(this == &v)
            return *this;
        else {
            x = v.x;
            y = v.y;
            A = v.A;
            v.x = v.y = 0;
            v.A = {};
        }
    }
    Matrix::~Matrix(){
        for (auto i : A)
        {
            i.clear();
        }
        A.clear();
        x = y = 0;
    }
    
    Matrix & Matrix::operator+=(const Matrix &v){
        if(x!=v.x || y!=v.y){
            throw std::runtime_error("+= : Matrixs are not Compatible");
        }else {
            for (auto i = 0; i < x; ++i){
                for (auto j = 0; j < y; ++j){
                    A[i][j] += v[i][j];
                }
            }
        }
        return *this;
    }
    Matrix & Matrix::operator-=(const Matrix &v){
        if(x!=v.x || y!=v.y){
            throw std::runtime_error("-= : Matrixs are not Compatible");
        }else {
            for (int i = 0; i < x; ++i)
            {
                for (int j = 0; j < y; ++j)
                {
                    A[i][j] -= v[i][j];
                }
            }
        }
        return *this;
    }
    
    Matrix operator*(const Matrix &v1, const Matrix &v2)
    {
        if(v1.y!=v2.x){
            throw std::runtime_error("* : Matrixs are not Compatible");
        }else {
            Matrix t(v1.x, v2.y);
            for (int i = 0; i < t.x; ++i){
                for (int j = 0; j < t.y; ++j){
                    t.A[i][j] = Util::innerProduct(v1[i], v2.columnAt(j));
                }
            }
            return t;
        }
    }
    Matrix operator+(const Matrix &v1, const Matrix &v2){
        Matrix t(v1);
        return t += v2;
    }
    Matrix operator-(const Matrix &v1, const Matrix &v2){
        Matrix t(v1);
        return t -= v2;
    }
    bool operator==(const Matrix &v1, const Matrix &v2){
        if(v1.x!=v2.x || v1.y!=v2.y)
            return false;
        for (int i = 0; i < v1.x; ++i){
            for (int j = 0; j < v1.y; ++j){
                if(v1[i][j]!=v2[i][j])
                    return false;
            }
        }
        return true;
    }
    bool operator!=(const Matrix &v1, const Matrix &v2){
        return !(v1 == v2);
    }

    Matrix Matrix::cheng(double a) const
    {
        Matrix t(*this);
        for (auto i = 0; i < x; ++i)
        {
            for (auto j = 0; j < y; ++j)
            {
                t[i][j] = t[i][j].cheng(a);
            }
        }
        return t;
    }
    Matrix Matrix::cheng(const Value &a) const{
        Matrix t(*this);
        for (auto i = 0; i < x; ++i)
        {
            for (auto j = 0; j < y; ++j)
            {
                t[i][j] *= a;
            }
        }
        return t;
    }
    std::string Matrix::toString()
    {
        std::stringstream st;
        st << " HAVE " << x << " rows, " << y << " columns:\n";
        for (long i = 0; i < x; ++i){
            st << "[ ";
            for (long j = 0; j < y; ++j){
                st << A[i][j].toString() << ", ";
            }
            st << "]\n";
        }
        return st.str();
    }
    unsigned int Matrix::rank() const {
        Matrix t(*this);
        // i-行 s-列 k-绝对值最大的行
        long i=0, k, s=0, j, m;
        vector<Value> v;
        Value e(0, 1);
        for (; i < x, s < y; ++i, ++s){
            k = i;
            for (; s < y; ++s){
                v = t.columnAt(s);
                k = Util::absmaxInVector(v ,i);
                if(v[k]!=Value(0,1))
                    break;
            } //s列之前的列均为0 s列的k行绝对值最大
            if(s==y) //i行i列右下方的元素全为0
                break;
            if(k > i) { //部分主元法 交换绝对值最大的行到前面
                swap(t.A[i], t.A[k]);
            }
            // 将i行以下 s列以后的元素设为新值
            for (j = i + 1; j < x; ++j){
                if(t[j][s]==Value(0,1)){
                    continue;
                }
                e = t[j][s] / t[i][s];
                for (m = s; m < y; ++m){
                    t[j][m] -= e * t[i][m];
                }
            }
        }
        return i;
    }

    // 初等变换矩阵 交换两行(列) 左乘该矩阵是做行变换
    Matrix Matrix::getE1(unsigned int size, unsigned int a, unsigned int b){
        Matrix t = Matrix::getI(size);
        swap(t[a], t[b]);
        return t;
    }
    // 初等变换矩阵 将某一行(列)乘以一个实数
    Matrix Matrix::getE2(unsigned int size, unsigned int a, double d)
    {
        return Matrix::getE3(size, a, a, d);
    }
    // 初等变换矩阵 将a行乘以一个实数加到b行上 a b应为不同的数
    Matrix Matrix::getE3(unsigned int size, unsigned int a, unsigned int b, double d)
    {
        Matrix t = Matrix::getI(size);
        t[b][a] = Value(d, 1);
        return t;
    }
    // LU分解 返回0分解成功 1不能对A分解 2分解失败
    int Matrix::LU(Matrix &P, Matrix &L, Matrix &U) const{
        if(!isSquare()){
            return 1;
        }
        P = getI(x);
        L = getI(x);
        U = Matrix(*this);
        // i-行 s-列 k-绝对值最大的行
        long i = 0, k, s = 0, j, m;
        vector<Value> v;
        Value e(0, 1);
        for (; i < x, s < x; ++i, ++s)
        {
            k = i;
            for (; s < x; ++s)
            {
                v = U.columnAt(s);
                k = Util::absmaxInVector(v, i);
                if (v[k] != Value(0, 1))
                    break;
            }           // s列之前的列均为0 s列的k行绝对值最大
            if (s == x) // i行i列右下方的元素全为0
                break;
            if (k > i)
            { //部分主元法 交换绝对值最大的行到前面
                swap(U.A[i], U.A[k]);
                swap(P.A[i], P.A[k]);
            }
            // 将i行以下 s列以后的元素设为新值
            for (j = i + 1; j < x; ++j)
            {
                if (U[j][s] == Value(0, 1))
                {
                    continue;
                }
                e = U[j][s] / U[i][s];
                U[j][s] = e;
                //std::cout << j << " , " << i << " : " << e.toString() << "\n";
                for (m = s+1; m < x; ++m)
                {
                    U[j][m] -= e * U[i][m];
                }
            }
            //std::cout << "U:" << U.toString();
        }
        if(i!=x)
            return 2;
        else {
            for (j = 0; j < x; ++j){
                for (m = 0; m < j; ++m){
                    L[j][m] = U[j][m];
                    U[j][m] = Value(0, 1);
                }
            }
            return 0;
        }
    }
    // QR分解 返回0分解成功 1不能对A分解 2分解失败
    int Matrix::QR(Matrix &Q, Matrix &R) const{
        if(rank()!=y){
            return 1;
        }
        Q = *this;
        R = getI(y);
        long k = 0,i,j;
        vector<Value> e,s;
        for (; k < y; ++k){
            i = 0;
            for (; i < k; ++i){
                R[i][k] = Util::innerProduct(Q.columnAt(i), columnAt(k));
            }
            j = 0;
            e = columnAt(k);
            for (; j < k; ++j){ //x(k)减去一系列投影 得到u(k)
                s = Q.columnAt(j);
                Util::chengVector(s, R[j][k]);
                for (int i = 0; i < e.size(); ++i){
                    e[i] -= s[i];
                }
            }
            R[k][k] = Util::innerProduct(e, e).sqrt();
            for (int t = 0; t < x; ++t){
                Q[t][k] = e[t] / R[k][k];
            }
        }
        return 0;
    }
    // Householder约简 返回0分解成功
    int Matrix::Householder(Matrix &P, Matrix &T) const{
        int n;
        if(y >= x)
            n = x - 1;
        else
            n = y;
        T = *this;
        P = getI(x);
        Matrix RR, R, I;
        int i = 0;
        vector<Value> u;
        for (; i < n; ++i){
            u = T.columnAt(i,i);
            I = getI(u.size());
            u[0] = u[0] - Util::innerProduct(u, u).sqrt();
            if(!Util::isAllZero(u)){
                R = I - (Matrix(u, 0) * Matrix(u, 1))
                            .cheng(Value(1, 1) / Util::innerProduct(u, u))
                            .cheng(2);
            }else {
                R = I;
            }
            RR = getI(x);
            for (int k = i; k < x; ++k){
                for (int j = i; j < x; ++j){
                    RR[k][j] = R[k - i][j - i];
                }
            }
            T = RR * T;
            P = RR * P;
        }
        return 0;
    }
    // Givens约简 返回0分解成功
    int Matrix::Givens(Matrix &P, Matrix &T) const{
        T = *this;
        P = getI(x);
        Matrix PP;
        Value c, s ,d;
        vector<Value> t;
        int i = 0, j = 0;
        for (; i < y-1; ++i){
            for (j=i+1; j < x; ++j){
                t = T.columnAt(i);
                d = (t[i].square() + t[j].square()).sqrt();
                c = t[i] / Value(d.toDouble(), 1);
                s = t[j] / Value(d.toDouble(), 1);
                PP = getI(x);
                PP[i][i] = PP[j][j] = c;
                PP[i][j] = s;
                PP[j][i] = s.cheng(-1);
                P = PP * P;
                T = PP * T;
            }
        }
        return 0;
    }
    // URV分解 返回0分解成功
    int Matrix::URV(Matrix &U, Matrix &R, Matrix &V) const{
        Matrix P1, T1, P2, T2;
        int r = rank();
        Givens(P1, T1);
        Matrix AT = T();
        AT.Givens(P2, T2); //使用Householder约简可能计算溢出
        U = Matrix(x, x);
        V = Matrix(y, y);
        Matrix V1 = fun1(T1); //V的前r列
        Matrix U1 = fun1(T2); //U的前r列
        if(V1.y != r || U1.y != r){
             throw std::runtime_error("V1 or U1 donot have r vectors");
        }
        for (int i = 0; i < r; ++i){
            for (int m = 0; m < x; ++m){
                U[m][i] = U1[m][i];
            }
            for (int n = 0; n < y; ++n){
                V[n][i] = V1[n][i];
            }
        }
        for (int i = r; i < x; ++i){
            for (int m = 0; m < x; ++m){
                U[m][i] = P1[i][m];
            }
        }
        for (int i = r; i < y;++i){
            for (int m = 0; m < y; ++m){
                V[m][i] = P2[i][m];
            }
        }
        R = U.T() * (*this) * V;
        return 0;
    }

    // T为行阶梯阵 将其中不为0的行作为一组基 然后转为标准正交基 按列存储
    Matrix Matrix::fun1(const Matrix &T) const{
        Matrix t(T.y, 0);
        vector<Value> v;
        for (int i = T.x - 1; i >= 0; --i){
            v = (T.A)[i];
            if(!Util::isAllZero(v)){
                for (int j = 0; j < T.y; ++j){
                    t[j].push_back(v[j]);
                }
                t.y++;
            }
        }
        Matrix Q, R;
        t.QR(Q, R);
        return Q;
    }
}
