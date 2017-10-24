//
// Created by gritukan on 10/24/17.
//

#ifndef GAUSS_REALNUMBER_H
#define GAUSS_REALNUMBER_H

#include <cmath>
#include "FieldNumber.h"

template<typename T>
class RealNumber : public FieldNumber {
private:
    T value;
public:
    RealNumber() = default;

    RealNumber(T x) {
        value = x;
    }

    operator T() const {
        return value;
    }

    RealNumber operator+(const RealNumber &other) {
        return value + other.value;
    }

    RealNumber operator*(const RealNumber &other) {
        return value * other.value;
    }

    RealNumber operator<(const RealNumber &other) {
        std::isless(value, other.value);
    }

    RealNumber operator==(const RealNumber &other) {
        return !std::isless(value, other.value) && !std::isless(other.value, value);
    }

    template <typename U>
    friend RealNumber<U> negate(const RealNumber<U> &x);

    template <typename U>
    friend RealNumber<U> inverse(const RealNumber<U> &x);

    static RealNumber<T> ONE, ZERO;
};

template <typename T>
RealNumber<T> RealNumber<T>::ONE = 1.0;

template <typename T>
RealNumber<T> RealNumber<T>::ZERO = 0.0;

template <typename T>
RealNumber<T> negate(const RealNumber<T> &x) {
    return RealNumber<T>(-x.value);
}

template <typename T>
RealNumber<T> inverse(const RealNumber<T> &x) {
    return RealNumber<T>(1.0 / x.value);
}

#endif //GAUSS_REALNUMBER_H
