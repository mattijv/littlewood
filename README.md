# Littlewood conjecture

This repo contains a program for computational investigation of the
[Littlewood conjecture](https://en.wikipedia.org/wiki/Littlewood_conjecture)
using continued fractions.


The code depends on [GMP](https://gmplib.org/) for the arbitrary precision
integer math, and on [Boost](https://www.boost.org/) for wide fixed width
integer math.

The code uses [{fmt}](https://github.com/fmtlib/fmt) for string formatting.

Check the makefile for options on compiling and running the program.

## main.cpp

Entrypoint of the program. Parses the CLI options, creates initial pairs
based on the selected option, initializes and executes the configured
number of threads to run the calculation. Prints out timing information.

## fractions.hpp

Contains types and helper methods related to the fractional types
(rational numbers, convergents and convergent pairs) used in the code,
as well as the logic to create the initial list of convergent pairs.

## littlewood.hpp

Contains the main part of the algorithm: the code to check if a
convergent pair meets the Littlewood criteria.

## modular_math.hpp

Contains helper functions to perform "modular" math.

## visualisation.hpp

Contains a helper method to print a textual representation of a
convergent pair. Used mostly for debugging.
