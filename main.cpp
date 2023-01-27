#include <iostream>
#include <span>
#include <thread>
#include <random>
#include <string>
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include "fractions.hpp"
#include "littlewood.hpp"

#include <boost/multiprecision/cpp_int.hpp>
using BigInt = boost::multiprecision::int1024_t;

template <typename T>
void check_pair(const fractions::convergent_pair<T>& pair, int depth, int N) {
    if (LW::meets_littlewood_criteria(pair, N)) {
        return;
    }

    // Napkin math assumption of a safe limit where we won't overflow the 1024bit fixed width int
    assert(depth < 80);

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
};

configuration parse_cli_arguments(int argc, char* argv[]) {
    configuration config = {10, 1};
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            auto argument = std::string(argv[i]);
            if (argument.substr(0,2) == "-N") {
                config.N = std::stoi(argument.substr(2));
            } else if (argument.substr(0,2) == "-j") {
                config.n_threads = std::stoi(argument.substr(2));
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

    auto pairs = fractions::convergent_pairs<BigInt>(config.N, 3);
    std::cout << fmt::format("Initial pairs: {}", pairs.size()) << std::endl;

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