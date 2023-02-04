#include <iostream>
#include <span>
#include <thread>
#include <random>
#include <string>
#include <cassert>
#include <cmath>
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include "fractions.hpp"
#include "littlewood.hpp"

#ifdef FIXED_WIDTH_INTEGERS
#ifndef INTEGER_WIDTH
#define INTEGER_WIDTH 1024
#endif
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/tools/precision.hpp>
#if INTEGER_WIDTH == 128
using BigInt = boost::multiprecision::int128_t;
#elif INTEGER_WIDTH == 256
using BigInt = boost::multiprecision::int256_t;
#elif INTEGER_WIDTH == 512
using BigInt = boost::multiprecision::int512_t;
#elif INTEGER_WIDTH == 1024
using BigInt = boost::multiprecision::int1024_t;
#endif
#endif

#ifndef FIXED_WIDTH_INTEGERS
#include <NTL/ZZ.h>
#define ARBITRARY_WIDTH_INTEGERS
// Use arbitrary integers instead.
using BigInt = NTL::ZZ;
#endif


template <typename T>
void check_pair(const fractions::convergent_pair<T>& pair, int depth, int N) {
    if (LW::meets_littlewood_criteria(pair, N)) {
        return;
    }

    #if defined(FIXED_WIDTH_INTEGERS) && defined(OVERFLOW_PROTECTION)
    // If we are using fixed width integers, this check should guarantee that we can't overflow the
    // integer in the next iteration.
    assert(boost::multiprecision::pow(pair.alpha.current.den, 5) < boost::math::tools::max_value<BigInt>() / (static_cast<BigInt>(2 * std::pow(N, 8))));
    #endif

    std::vector<fractions::convergent_pair<T>> child_pairs = {};
    fractions::subdivide(pair, N, child_pairs);
    for(auto& child_pair: child_pairs) {
        check_pair(child_pair, depth + 1, N);
    }
}

template <typename T>
void check_pair_list(const std::span<fractions::convergent_pair<T>>& pairs, int N) {
    for (auto& pair: pairs) {
        check_pair(pair, 0, N);
    }
}

struct configuration {
    int N;
    uint n_threads;
    uint subdivisions;
};

configuration parse_cli_arguments(int argc, char* argv[]) {
    configuration config = {10, 1, 0};
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            auto argument = std::string(argv[i]);
            if (argument.substr(0,2) == "-N") {
                config.N = std::stoi(argument.substr(2));
            } else if (argument.substr(0,2) == "-j") {
                config.n_threads = std::stoi(argument.substr(2));
            } else if (argument.substr(0,2) == "-s") {
                config.subdivisions = std::stoi(argument.substr(2));
            }
        }
    }

    assert(config.n_threads <= std::thread::hardware_concurrency());
    assert(config.N > 0);
    return config;
}

int main(int argc, char* argv[]) {

    auto config = parse_cli_arguments(argc, argv);

    std::cout << fmt::format("Running on {} thread(s) with N={}.", config.n_threads, config.N) << std::endl;

    #ifdef FIXED_WIDTH_INTEGERS
    std::cout << fmt::format("Using fixed width integers, with bit width of {}.", INTEGER_WIDTH) << std::endl;
    #endif

    #ifdef ARBITRARY_WIDTH_INTEGERS
    std::cout << fmt::format("Using arbitrary sized integers.") << std::endl;
    #endif

    auto pairs = fractions::convergent_pairs<BigInt>(config.N, config.subdivisions);
    std::cout << fmt::format("Initial pairs (from {} subdivisions): {}", config.subdivisions, pairs.size()) << std::endl;

    // Shuffle the pairs so that the different sized numerators are evenly split among the
    // threads. Otherwise one thread might end up processing all the pathological ones, and we
    // lose a lot of the value from multithreading.
    auto rd = std::random_device {};
    auto rng = std::default_random_engine {rd()};
    std::shuffle(std::begin(pairs), std::end(pairs), rng);

    std::vector<std::thread> threads;

    int pairs_per_thread = pairs.size() / config.n_threads;

    // If the number of pairs is not neatly divisible by the number of threads
    // the last thread will get the extra pairs, as calculated here.
    int rest = pairs.size() - (config.n_threads - 1) * pairs_per_thread;
    
    for (uint i = 0; i < config.n_threads; i++) {
        int start = i * pairs_per_thread;
        int number_of_pairs_to_process = i < config.n_threads ? pairs_per_thread : rest;
        threads.emplace_back(
            check_pair_list<BigInt>,
            std::span<fractions::convergent_pair<BigInt>>{pairs}.subspan(start, number_of_pairs_to_process),
            config.N
        );
    }

    for (std::size_t i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    std::cout << "Done" << std::endl;

    return 0;
}
