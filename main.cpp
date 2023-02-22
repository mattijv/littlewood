#include <iostream>
#include <span>
#include <thread>
#include <random>
#include <string>
#include <cassert>
#include <cmath>
#include <mutex>
#include <functional>
#include <chrono>
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

std::mutex work_queue_mutex;
std::mutex thread_status_mutex;

int thread_count = 1;
int done_thread_count = 0;

struct configuration {
    int N;
    uint n_threads;
    uint subdivisions;
    uint buckets;
    uint bucket;
};

// Very naive CLI argument parser
configuration parse_cli_arguments(int argc, char* argv[]) {
    configuration config = {10, 1, 0, 1, 1};
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            auto argument = std::string(argv[i]);
            auto handle = argument.substr(0,2);
            if (handle == "-N") {
                config.N = std::stoi(argument.substr(2));
            } else if (handle == "-j") {
                config.n_threads = std::stoi(argument.substr(2));
            } else if (handle == "-s") {
                config.subdivisions = std::stoi(argument.substr(2));
            } else if (handle == "-B") {
                config.buckets = std::stoi(argument.substr(2));
            } else if (handle == "-b") {
                config.bucket = std::stoi(argument.substr(2));
            }
        }
    }

    assert(config.n_threads <= std::thread::hardware_concurrency());
    assert(config.N > 0);
    assert(config.buckets >= 1);
    assert(config.bucket > 0 && config.bucket <= config.buckets);
    return config;
}

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

bool all_threads_done() {
    thread_status_mutex.lock();
    bool all_done = thread_count == done_thread_count;
    thread_status_mutex.unlock();
    return all_done;
}

void increment_done_thread_count() {
    thread_status_mutex.lock();
    done_thread_count++;
    thread_status_mutex.unlock();
}

void decrement_done_thread_count() {
    thread_status_mutex.lock();
    done_thread_count--;
    thread_status_mutex.unlock();
}

template <typename T>
void process(std::vector<fractions::convergent_pair<T>>& queue, int N) {
    
    bool done = false;

    // Until all threads have marked themselves as done, we need to try processing pairs from the queue.
    while (!all_threads_done()) {
        
        work_queue_mutex.lock();
        // If there is work available
        if (queue.size() > 0) {
            // If we previously marked ourself to be done:
            if (done) {
                // Decrease the count of done threads by one
                decrement_done_thread_count();
                // and mark ourself as not done.
                done = false;
            }

            // Get a new pair to work on and unlock the mutex.
            auto pair = std::move(queue.back());
            queue.pop_back();
            work_queue_mutex.unlock();

            // The pair passes the criteria, so we can forget about it and jump back to the start of the loop.
            if (LW::meets_littlewood_criteria(pair, N)) {
                continue;
            }

            #if defined(FIXED_WIDTH_INTEGERS) && defined(OVERFLOW_PROTECTION)
            // If we are using fixed width integers, this check should guarantee that we can't overflow the
            // integer in the next iteration.
            assert(boost::multiprecision::pow(pair.alpha.current.den, 5) < boost::math::tools::max_value<BigInt>() / (static_cast<BigInt>(2 * std::pow(N, 8))));
            #endif

            // The pair does not match the criteria, so we divide it into N new pairs with larger denominators
            std::vector<fractions::convergent_pair<T>> child_pairs = {};
            fractions::subdivide(pair, N, child_pairs);
            // and add them to the work queue.
            work_queue_mutex.lock();
            queue.insert(std::end(queue), std::begin(child_pairs), std::end(child_pairs));
            work_queue_mutex.unlock();
        } else {
            // No work available on the queue
            work_queue_mutex.unlock();
            // If we have not yet marked ourself as done:
            if (!done) {
                // Mark ourself as done
                done = true;
                // and increment the number of done threads.
                increment_done_thread_count();
            }
            // Sleep for 2 seconds and before continuing the while loop.
            using namespace std::chrono_literals;
            std::this_thread::sleep_for(2000ms);
        }
    }
}

template <typename T>
std::vector<fractions::convergent_pair<T>> select_bucket(std::vector<fractions::convergent_pair<T>>& pairs, configuration config) {
    assert(pairs.size() >= config.buckets);

    if (config.buckets == 1) {
        return pairs;
    }

    std::vector<fractions::convergent_pair<T>> bucket = {};

    int index = config.bucket - 1;
    while (index < pairs.size()) {
        bucket.push_back(pairs[index]);
        index += config.buckets;
    }

    return bucket;
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

    auto bucket = select_bucket<BigInt>(pairs, config);
    std::cout << fmt::format("Pairs split into {} buckets, processing bucket #{} which has {} pairs.", config.buckets, config.bucket, bucket.size()) << std::endl;

    thread_count = config.n_threads;

    std::vector<std::thread> threads;
    for (uint i = 0; i < config.n_threads; i++) {
            threads.emplace_back(
                process<BigInt>,
                std::ref(bucket),
                config.N
            );
    }

    for (std::size_t i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    std::cout << "Done" << std::endl;

    return 0;
}
