#pragma once
#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <string>
#include <vector>
#include <sstream>
#include "base.h"

using std::abs;
using std::string;
using std::swap;
using std::vector;

class MatrixWork::Util{
    private:
        Util();
    public:
        static const double error;
        static int gcd(int a, int b);
        static bool isInteger(double a){
            if(a-long(a) < Util::error)
                return true;
            else
                return false;
        }
        static bool isAllZero(const vector<Value> &v);
        static string toS(double a){
            std::stringstream ss;
            ss << a;
            return ss.str();
        }
        static Value innerProduct(const vector<Value>& v1, const vector<Value>& v2);
        static void chengVector(vector<Value> &v, const Value &a);
        static vector<Value> chengVector(unsigned int size, unsigned int n, const Value &a);
        static unsigned int absmaxInVector(const vector<Value> &v, int from = 0);
        static void testLU();
        static void testQR();
        static void testHouseholder();
        static void testGivens();
        static void testURV();
        static Value getValue(const string &v);
};

#endif