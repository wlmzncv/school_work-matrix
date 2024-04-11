#include <iostream>
#include "Base.h"
#include "Matrix.h"
#include "Solution.h"
#include <sstream>


using MatrixWork::Matrix;
using MatrixWork::Solution;
using MatrixWork::Util;
using MatrixWork::Value;
using std::cin;
using std::cout;
using std::endl;
using std::string;

static string tip = "============****************============\n";
static string tip1 = "Input 1 to set Matrix A\n";
static string tip2 = "Input 2 to set b of Linear Equations\n";
static string tip3 = "Input 3 to take LU Factorization\n";
static string tip4 = "Input 4 to take QR Factorization\n";
static string tip5 = "Input 5 to take Householder Reduction\n";
static string tip6 = "Input 6 to take Givens Reduction\n";
static string tip7 = "Input 7 to take URV Factorization\n";
static string tip11 = "Input 8 to get A's Determinant\n";
static string tip0 = "Input 0 to exit\n";
static string tip8 = "Input m(rows' number) and n(columns' number) for Matrix A\n";
static string tip9 = "Input elements for b\n";
static string tip10 = "Please set Matrix A first\n";

static Solution solution;

// 输入矩阵 无先决条件
void fun1(){
    cout << tip8;
    unsigned int m, n;
    cout << "m=";
    cin >> m;
    cout << "n=";
    cin >> n;
    solution.A = Matrix(m, n);
    Value v;
    string s;
    for (int i = 0; i < m; ++i){
        for (int j = 0; j < n; ++j){
            cout << "A[" << i + 1 << "][" << j + 1 << "]=";
            cin >> s;
            solution.A[i][j] = Util::getValue(s);
        }
    }
    if(solution.b.X() != m)
        solution.b = Matrix(m, 1);
    solution.freshState();
    cout << "A "
         << solution.A.toString();
}
// 输入b 需要先输入矩阵
void fun2(){
    if(solution.A.empty()){
        cout << tip10;
        return;
    }
    string s;
    cout << tip9;
    for (int i = 0; i < solution.A.X(); ++i){
        cout << "b[" << i + 1 << "]=";
        cin >> s;
        solution.b[i][0] = Util::getValue(s);
    }
    cout << "b "
         << solution.b.toString();
}
// 执行完一个操作后停顿 输入回车键后继续
void fun3(){
    string s;
    getline(cin, s);
    cout << "####### Please press Enter to Continue...";
    getline(cin, s);
}


int main(){
    unsigned int oper = -1;
    bool going = true;
    string rubbish;
    while(going){
        switch (oper)
        {
        case 0:
            going = false;
            cout << "Good Bye!" << endl;
            break;
        case 1:
            fun1();
            fun3();
            break;
        case 2:
            fun2();
            fun3();
            break;
        case 3:
            solution.takeLU();
            fun3();
            break;
        case 4:
            solution.takeQR();
            fun3();
            break;
        case 5:
            solution.takeHouseholder();
            fun3();
            break;
        case 6:
            solution.takeGivens();
            fun3();
            break;
        case 7:
            solution.takeURV();
            fun3();
            break;
        case 8:
            solution.getDeterminant();
            fun3();
            break;
        case 4294967295:
            break;
        default:
            cout << "wrong letter: " << oper << endl;
            fun3();
            break;
        }
        if(going){
            cout << tip;
            cout << "A " << solution.A.toString() << "b " << solution.b.toString() << solution.getState();
            cout << tip1 << tip2 << tip3 << tip4 << tip5 << tip6 << tip7 << tip11 << tip0;
            cout << "Choose:";
            cin >> oper;
            getline(cin, rubbish);
            cout << "You have choose " << oper << endl;
        }
    }
    return 0;
}