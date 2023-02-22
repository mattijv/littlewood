# Compilation options
compiler = g++
common_flags = -g -O2 -std=c++20 -pthread -march=native
gmp_libs = -lntl -lgmp
target = main.cpp
output = -o build/lw
bits ?= 256

# Run options
N ?= 9
threads ?= 1
subdivisions = 0
buckets ?= 1
bucket ?= 1

compile:
	$(compiler) $(common_flags) $(target) $(output) $(gmp_libs)

compile-fixed:
	$(compiler) $(common_flags) $(target) $(output) -DFIXED_WIDTH_INTEGERS -DINTEGER_WIDTH=$(bits)

compile-fixed-safe:
	$(compiler) $(common_flags) $(target) $(output) -DFIXED_WIDTH_INTEGERS -DINTEGER_WIDTH=$(bits) -DOVERFLOW_PROTECTION

_run:
	/usr/bin/time --format="Executed in %E" ./build/lw -N$(N) -j$(threads) -s$(subdivisions)

run: compile _run

run-fixed: compile-fixed _run

run-fixed-safe: compile-fixed-safe _run
