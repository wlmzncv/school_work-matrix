#include "Value.h"

namespace MatrixWork{
    Value& Value::operator=(const Value &v)
    {
        numerator = v.numerator;
        denominator = v.denominator;
        return *this;
    }
    Value& Value::operator+=(const Value &v)
    {
        numerator = numerator * v.denominator + v.numerator * denominator;
        denominator *= v.denominator;
        selfReduction();
        return *this;
    }
    Value operator+(const Value &v1, const Value &v2)
    {
        Value t(v1);
        t += v2;
        return t;
    }
    Value &Value::operator-=(const Value &v)
    {
        numerator = numerator * v.denominator - v.numerator * denominator;
        denominator *= v.denominator;
        selfReduction();
        return *this;
    }
    Value operator-(const Value &v1, const Value &v2)
    {
        Value t(v1);
        t -= v2;
        return t;
    }
    Value &Value::operator*=(const Value &v)
    {
        denominator *= v.denominator;
        numerator *= v.numerator;
        selfReduction();
        return *this;
    }
    Value operator*(const Value &v1, const Value &v2)
    {
        Value t(v1);
        t *= v2;
        return t;
    }
    Value &Value::operator/=(const Value &v)
    {
        denominator *= v.numerator;
        numerator *= v.denominator;
        selfReduction();
        return *this;
    }
    Value operator/(const Value &v1, const Value &v2)
    {
        Value t(v1);
        t /= v2;
        return t;
    }
    bool operator>(const Value &v1, const Value &v2)
    {
        return v1.toDouble() > v2.toDouble();
    }
    bool operator<(const Value &v1, const Value &v2)
    {
        return v1.toDouble() < v2.toDouble();
    }
    bool operator==(const Value &v1, const Value &v2)
    {
        return v1.toDouble() == v2.toDouble();
    }
    bool operator!=(const Value &v1, const Value &v2)
    {
        return v1.toDouble() != v2.toDouble();
    }

}

