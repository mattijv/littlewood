# Copyright 2023 Matti Vapa

# Dependencies
deps_directory = dependencies
fmt_path = $(deps_directory)/fmt
fmt_repo = git@github.com:fmtlib/fmt.git
boost_url = https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
boost_directory = boost_1_81_0
boost_archive = $(boost_directory).tar.gz
boost_path = $(deps_directory)/$(boost_directory)

# CSC variables
project = project_2007026
run_directory = /projappl/$(project)/littlewood

# Compilation options
compiler = g++
debug = -g
common_flags = -O3 -std=c++20 -pthread -march=native
gmp_libs = -lntl -lgmp
boost_lib = -I $(boost_path)
target = main.cpp
output = -o build/lw
bits ?= 256

# Run options
N ?= 9
threads ?= 1
subdivisions = 0
buckets ?= 1
bucket ?= 1

# Download dependencies
$(fmt_path):
	@mkdir -p $(deps_directory)
	git clone $(fmt_repo) $(fmt_path)

$(boost_path):
	@mkdir -p $(deps_directory)
	wget -P $(deps_directory)/ $(boost_url)
	cd $(deps_directory) && tar -xzf $(boost_archive) && rm $(boost_archive)

# Compilation targets
compile: $(fmt_path)
	$(compiler) $(debug) $(common_flags) $(target) $(output) $(gmp_libs)

compile-fixed: $(fmt_path) $(boost_path)
	$(compiler) $(debug) $(common_flags) $(boost_lib) $(target) $(output) -DFIXED_WIDTH_INTEGERS -DINTEGER_WIDTH=$(bits)

compile-fixed-safe: $(fmt_path) $(boost_path)
	$(compiler) $(debug) $(common_flags) $(boost_lib) $(target) $(output) -DFIXED_WIDTH_INTEGERS -DINTEGER_WIDTH=$(bits) -DOVERFLOW_PROTECTION

compile-production: $(fmt_path) $(boost_path)
	$(compiler) $(common_flags) $(boost_lib) $(target) $(output) -DFIXED_WIDTH_INTEGERS -DINTEGER_WIDTH=$(bits) -DOVERFLOW_PROTECTION
	mkdir -p $(run_directory)
	cp build/lw $(run_directory)/

# Run targets

_run:
	./build/lw -N$(N) -j$(threads) -s$(subdivisions) -B$(buckets) -b$(bucket)

run: compile _run

run-fixed: compile-fixed _run

run-fixed-safe: compile-fixed-safe _run

print-pairs: compile-fixed
	./build/lw -N$(N) -s$(subdivisions) -B$(buckets) -b$(bucket) -p
