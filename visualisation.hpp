/**
 * Copyright 2023 Topi Törmä, Matti Vapa
 */

#ifndef VISUALISATION
#define VISUALISATION

#include <sstream>
#include "fractions.hpp"

namespace visualisation {

    /**
     * This method will print a convergent pair in the format "{a/b, c/d}, {e/f, g/h}"
     */
    template <typename T>
    std::stringstream string_representation(const fractions::convergent_pair<T>& pair) {
        std::stringstream output;
        output << "{";
        output << pair.alpha.current.num << "/" << pair.alpha.current.den;
        output << ", ";
        output << pair.alpha.previous.num << "/" << pair.alpha.previous.den;
        output << "}, {";
        output << pair.beta.current.num << "/" << pair.beta.current.den;
        output << ", ";
        output << pair.beta.previous.num << "/" << pair.beta.previous.den;
        output << "}";
        return output;
    }
}

#endif