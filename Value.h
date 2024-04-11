#pragma once
#ifndef VALUE_H
#define VALUE_H

#include <stdexcept>
#include <string>
#include <sstream>
#include "Base.h"
#include "math.h"
#include "Util.h"

class MatrixWork::Value
{

private:
    static const int MAX_SIZE = 20;
    double numerator = 0; //分子
    double denominator = 1; //分母
    void selfReduction(){
        if(Util::isInteger(numerator)
        && Util::isInteger(denominator)){
            numerator = int(numerator);
            denominator == int(denominator);
            double r = Util::gcd(int(numerator), int(denominator));
            numerator /= r;
            denominator /= r;
        }
        if (numerator==0)
        {
            denominator = 1;
        }
        if(denominator < 0){
            denominator = -denominator;
            numerator = -numerator;
        }
    }

public:
    Value(){};
    Value(double a, double b);
    Value(const Value &v) : numerator(v.numerator), denominator(v.denominator) {}
    Value(double a):numerator(a),denominator(1){}
    Value(int a):numerator(a),denominator(1){}
    Value& operator=(const Value &v);

    Value &operator+=(const Value &v);
    Value &operator-=(const Value &v);
    Value &operator*=(const Value &v);
    Value &operator/=(const Value &v);

    bool isPositive() const
    {
        return numerator * denominator > 0;
    }
    Value abs() const {
        Value t(*this);
        if(numerator<0)
            t.numerator = -numerator;
        if(denominator<0)
            t.denominator = -denominator;
        return t;
    }
    Value square() const {
        Value t;
        t.numerator = numerator * numerator;
        t.denominator = denominator * denominator;
        t.selfReduction();
        return t;
    }
    Value cheng(double b) const {
        Value t(*this);
        t.numerator *= b;
        t.selfReduction();
        return t;
    }
    Value sqrt() const {
        Value t(*this);
        t.numerator = std::sqrt(numerator);
        t.denominator = std::sqrt(denominator);
        t.selfReduction();
        return t;
    }
    double toDouble() const {
        return numerator / denominator;
    }
    std::string toString() const{
        if(denominator == 1){
            return Util::toS(numerator);
        }else {
            std::string s = Util::toS(numerator) + "/" + Util::toS(denominator);
            if(s.length()>MAX_SIZE){
                return Value(toDouble(), 1).toString();
            }else
                return s;
        }
    }
    ~Value(){};
};

inline MatrixWork::Value::Value(double a, double b) :numerator(a), denominator(b)
{
    if(b==0){
        throw std::runtime_error("denominator of Fraction cannot be 0");
    }
}

#endif