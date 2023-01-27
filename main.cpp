#include <iostream>
#include <span>
#include <thread>
#include <random>
#include "fractions.hpp"
#include "littlewood.hpp"

#include <boost/multiprecision/cpp_int.hpp>
namespace mp = boost::multiprecision;
using BigInt = mp::int1024_t;

template <typename T>
void check_pair(const fractions::convergent_pair<T>& pair, int depth, int N) {
    if (LW::meets_littlewood_criteria(pair, N)) {
        return;
    }

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

int main(int argc, char* argv[]) {
    int N = 9;
    auto pairs = fractions::convergent_pairs<BigInt>(N, 3);
    std::cout << "Initial pairs: " << pairs.size() << std::endl;

    auto rd = std::random_device {};
    auto rng = std::default_random_engine {rd()};
    std::shuffle(std::begin(pairs), std::end(pairs), rng);


    int n_threads = 8;
    std::vector<std::thread> threads;

    int pairs_per_thread = pairs.size() / n_threads;

    // If the number of pairs is not neatly divisible by the number of threads
    // the last thread will get the extra pairs, calculated here.
    int rest = pairs.size() - (n_threads - 1) * pairs_per_thread;

    // Initialize all threads but the last.
    
    for (int i = 0; i < n_threads; i++) {
        int start = i * pairs_per_thread;
        int number_of_pairs_to_process = i < n_threads ? pairs_per_thread : rest;
        threads.emplace_back(check_pair_list<BigInt>, std::span<fractions::convergent_pair<BigInt>>{pairs}.subspan(start, number_of_pairs_to_process), N);
    }

    for (uint i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    std::cout << "Done" << std::endl;
}