/**
 * Copyright 2023 Topi Törmä, Matti Vapa
 */

#ifndef MODULAR_MATH
#define MODULAR_MATH

#include <algorithm>
#include "fractions.hpp"

namespace modular_math {
    template <typename T>
    T remainder_with_least_absolute_value(T q, fractions::rational<T> number) {
        T remainder = (q * number.num) % number.den;
        return std::min(remainder, number.den - remainder);
    }

    template <typename T>
    T modular_addition(T augend, T addend, T modulus) {
        T sum = augend + addend;
        if (sum <= modulus) {
            return sum;
        } else {
            return sum - modulus;
        }
    }
    
    template <typename T>
    T modular_substraction(T minuend, T subtrahend, T modulus) {
        T diff = minuend - subtrahend;
        if (diff >= 0) {
            return diff;
        } else {
            return diff + modulus;
        }
    }
}

#endif