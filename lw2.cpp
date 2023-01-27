#include <iostream>
#include <vector>
#include <span>
#include <thread>
#include <algorithm>
#include <random>
#include <NTL/ZZ.h>

#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;

//#define TRACK_CUTOFF_LOCATION
//#define PRINT_ENABLED
//#define TRACK_LARGEST
#define TRACK_MAX_DEPTH

#define BIGINT mp::int1024_t
#define i(number) static_cast<BIGINT>(number)
//#define i(number) NTL::ZZ(number)
//#define BIGINT NTL::ZZ
typedef BIGINT bigint;

long N = 9;

#ifdef TRACK_LARGEST
bigint chonker = i(0);
#endif

#ifdef TRACK_MAX_DEPTH
int max_depth = 0;
#endif

struct rational {
    bigint num, den;
};

struct convergent {
    rational current, previous;
};

struct convergent_pair {
    convergent alpha, beta;
};

convergent next_convergent(const convergent& base, bigint digit) {
    return {
        {
            digit * base.current.num + base.previous.num,
            digit * base.current.den + base.previous.den
        },
        base.current
    };
};

bigint littlewood(bigint Q, bigint a, bigint b) {
    bigint result = 2 * N * Q * (Q + a) * (Q + b);
    #ifdef TRACK_LARGEST
    chonker = std::max(chonker, result);
    #endif
    return result;
};

bigint remainder_with_least_absolute_value(bigint q, rational number) {
    bigint remainder = (q * number.num) % number.den;
    return std::min(remainder, number.den - remainder);
};

bigint modular_addition(bigint augend, bigint addend, bigint modulus) {
    bigint sum = augend + addend;
    if (sum <= modulus) {
        return sum;
    } else {
        return sum - modulus;
    }
};

bigint modular_substraction(bigint minuend, bigint subtrahend, bigint modulus) {
    bigint diff = minuend - subtrahend;
    if (diff >= 0) {
        return diff;
    } else {
        return diff + modulus;
    }
}

void debug_location(int location) {
    std::cout << location << std::endl;
}

#ifdef PRINT_ENABLED
void print_convergent(const convergent& conv) {
    std::cout << conv.current.num << "/" << conv.current.den << "\n";
    std::cout << conv.previous.num << "/" << conv.previous.den << "\n";
}
#endif

bool pair_meets_littlewood_criteria(const convergent_pair& pair) {
    
    convergent alpha = pair.alpha;
    convergent beta = pair.beta;

    bigint alpha_sum = alpha.current.den + alpha.previous.den;
    bigint beta_sum = beta.current.den + beta.previous.den;

    bigint epsilon = alpha.current.den * alpha_sum * beta.current.den * beta_sum;

    #ifdef TRACK_LARGEST
    chonker = std::max(chonker, epsilon);
    #endif
    
    //Quickly check if the bigger convergent denominator is good
    bigint alpha_remainder = remainder_with_least_absolute_value(beta.current.den, alpha.current);
    if (littlewood(beta.current.den, i(0), alpha_remainder * (alpha.current.den + alpha.previous.den)) < epsilon) {
        #ifdef TRACK_CUTOFF_LOCATION
        debug_location(1);
        #endif
        return true;
    }

    bigint alpha_diff = alpha.current.den - alpha.previous.den;
    bigint beta_diff = beta.current.den - beta.previous.den;

    bigint denominators[4] = {
        alpha.previous.den,
        alpha_diff,
        beta.previous.den,
        beta_diff
    };

    bigint alpha_den_beta_num = (alpha.previous.den * beta.current.num) % beta.current.den;
    bigint alpha_den_comp_beta_num = (alpha_diff * beta.current.num) % beta.current.den;
    bigint beta_den_alpha_num = (beta.previous.den * alpha.current.num) % alpha.current.den;
    bigint beta_den_comp_alpha_num = (beta_diff * alpha.current.num) % alpha.current.den;
    bigint multiple_alpha = (alpha.current.den * beta.current.num) % beta.current.den;

    bigint ab_rem = alpha_den_beta_num;
    bigint acb_rem = alpha_den_comp_beta_num;
    bigint ba_rem = beta_den_alpha_num;
    bigint bca_rem = beta_den_comp_alpha_num;
    bigint ma_rem = multiple_alpha;

    bigint max_remainder = 1 + std::max(i(N), pair.beta.previous.den / (N * N));
    bigint target_remainder = i(1);
    
    while (target_remainder < max_remainder) {

        if (littlewood(denominators[0], target_remainder * alpha_sum, std::min(ab_rem, beta.current.den - ab_rem) * beta_sum) < epsilon) {
            #ifdef TRACK_CUTOFF_LOCATION
            debug_location(2);
            #endif
            return true;
        }

        if (littlewood(denominators[1], target_remainder * alpha_sum, std::min(acb_rem, beta.current.den - acb_rem) * beta_sum) < epsilon) {
            #ifdef TRACK_CUTOFF_LOCATION
            debug_location(3);
            #endif
            return true;
        }

        if (littlewood(denominators[2], std::min(ba_rem, alpha.current.den - ba_rem) * alpha_sum, target_remainder * beta_sum) < epsilon) {
            #ifdef TRACK_CUTOFF_LOCATION
            debug_location(4);
            #endif
            return true;
        }

        if (littlewood(denominators[3], std::min(bca_rem, alpha.current.den - bca_rem) * alpha_sum, target_remainder * beta_sum) < epsilon) {
            #ifdef TRACK_CUTOFF_LOCATION
            debug_location(5);
            #endif
            return true;
        }

        if (littlewood(target_remainder * alpha.current.den, i(0), std::min(ma_rem, beta.current.den - ma_rem) * beta_sum) < epsilon) {
            #ifdef TRACK_CUTOFF_LOCATION
            debug_location(6);
            #endif
            return true;
        }

        target_remainder++;

        bigint ab_sum = denominators[0] + alpha.previous.den;

        if (ab_sum <= alpha.current.den) {
            denominators[0] = ab_sum;
            ab_rem = modular_addition(ab_rem, alpha_den_beta_num, beta.current.den);
        } else {
            denominators[0] = ab_sum - alpha.current.den;
            ab_rem = modular_substraction(ab_rem, alpha_den_comp_beta_num, beta.current.den);
        }

        bigint acb_sum = denominators[1] + alpha_diff;
        if (acb_sum <= alpha.current.den) {
            denominators[1] = acb_sum;
            acb_rem = modular_addition(acb_rem, alpha_den_comp_beta_num, beta.current.den);
        }
        else {
            denominators[1] = acb_sum - alpha.current.den;
            acb_rem = modular_substraction(acb_rem, alpha_den_beta_num, beta.current.den);
        }

        bigint ba_sum = denominators[2] + beta.previous.den;
        if (ba_sum <= beta.current.den) {
            denominators[2] = ba_sum;
            ba_rem = modular_addition(ba_rem, beta_den_alpha_num, alpha.current.den);
        }
        else {
            denominators[2] = ba_sum - beta.current.den;
            ba_rem = modular_substraction(ba_rem, beta_den_comp_alpha_num, alpha.current.den);
        }

        bigint bca_sum = denominators[3] + beta_diff;
        if (bca_sum <= beta.current.den) {
            denominators[3] = bca_sum;
            bca_rem = modular_addition(bca_rem, beta_den_comp_alpha_num, alpha.current.den);
        }
        else {
            denominators[3] = bca_sum - beta.current.den;
            bca_rem = modular_substraction(bca_rem, beta_den_alpha_num, alpha.current.den);
        }

        ma_rem = modular_addition(ma_rem, multiple_alpha, beta.current.den);
    }
    
    return false;
}

void subdivide_and_append_new_pairs(const convergent_pair& pair, std::vector<convergent_pair>& new_pairs) {
    convergent primary, secondary;
    if (pair.alpha.current.den < pair.beta.current.den) {
        primary = pair.alpha;
        secondary = pair.beta;
    } else {
        primary = pair.beta;
        secondary = pair.alpha;
    }

    for (int i = 1; i < N; i++) {
        convergent next = next_convergent(primary, i(i));
        
        if (next.current.den < secondary.current.den) {
            new_pairs.push_back({next, secondary});
        } else {
            new_pairs.push_back({secondary, next});
        }
    }
}

std::vector<convergent> initial_convergents() {
    std::vector<convergent> convergents = {};
    
    for (int i = 2; i < N; i++) {
        convergents.push_back({
            {i(1),i(i)},
            {i(0),i(1)}
        });
    }
    
    return convergents;
}

std::vector<convergent_pair> create_convergent_pairs() {
    std::vector<convergent> convergents = {};
    
    for (int i = 2; i < N; i++) {
        convergents.push_back({
            {i(1),i(i)},
            {i(0),i(1)}
        });
    }
    
    std::vector<convergent_pair> pairs = {};
    
    for (uint i = 0; i < convergents.size(); i++) {
        for (uint j = 0; j <= i; j++) {
            convergent a = convergents[i];
            convergent b = convergents[j];
            
            if (a.current.den * b.current.den >= 2 * N) {
                continue;
            }

            if (a.current.den < b.current.den) {
                pairs.push_back({a, b});
            } else {
                pairs.push_back({b, a});
            }
        }
    }

    for (int i = 0; i <= 2; i++) {
        std::vector<convergent_pair> next = {};
        for (auto& pair: pairs) {
            subdivide_and_append_new_pairs(pair, next);
        }
        swap(next, pairs);
    }
    
    return pairs;
}

void check_pair(const convergent_pair& pair, int depth) {
    if (pair_meets_littlewood_criteria(pair)) {
        #ifdef TRACK_MAX_DEPTH
        max_depth = std::max(max_depth, depth);
        #endif
        return;
    }

    assert(depth < 80);

    // get best Q from pair_meets_littlewood_criteria
    // in subdivide_etc calculate for each i > 1
    // new epsilon (add (i - 1) * alpha.current.den to alpha_sum, calculate epsilon as before)
    // A = remainder_with_least_absolute_value(bestQ * primary.current.num, primary.current.den)
    // B = remainder_with_least_absolute_value(bestQ * secondary.current.num, secondary.current.den)
    // if littlewood(bestQ, A * new_alpha_sum, B * beta_sum) < new_epsilon, no need for child or any of the following

    std::vector<convergent_pair> child_pairs = {};
    //child_pairs.reserve(N);
    subdivide_and_append_new_pairs(pair, child_pairs);
    for(auto& child_pair: child_pairs) {
        check_pair(child_pair, depth + 1);
    }
}

void check_pair_list(const std::span<convergent_pair>& pairs) {
    for (auto& pair: pairs) {
        check_pair(pair, 0);
    }
}

void echo(std::string text) {
    std::cout << text << std::endl;
}

int main(int argc, char* argv[]) {
    
    //po::options_description desc("Options");
    //desc.add_options()("threads,t", po::value<int>(), "number of threads to use");
    /*
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    int n_threads = vm.count("threads") ? vm["threads"].as<int>() : 1;
    */
    int n_threads = 8;
    auto pairs = create_convergent_pairs();

    auto rd = std::random_device {};
    auto rng = std::default_random_engine {rd()};
    std::shuffle(std::begin(pairs), std::end(pairs), rng);
    #ifndef TRACK_CUTOFF_LOCATION
    //int n_pairs = pairs.size();
    std::cout << "Initial pairs: " << pairs.size() << std::endl;
    #endif

    std::vector<std::thread> threads;

    int pairs_per_thread = pairs.size() / n_threads;

    // If the number of pairs is not neatly divisible by the number of threads
    // the last thread will get the extra pairs, calculated here.
    int rest = pairs.size() - (n_threads - 1) * pairs_per_thread;

    // Initialize all threads but the last.
    
    for (int i = 0; i < n_threads; i++) {
        int start = i * pairs_per_thread;
        int number_of_pairs_to_process = i < n_threads ? pairs_per_thread : rest;
        threads.emplace_back(check_pair_list, std::span<convergent_pair>{pairs}.subspan(start, number_of_pairs_to_process));
    }

    for (uint i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    /*
    int handled_pairs = 0;
    for (auto& pair: pairs) {
        check_pair(pair);
        handled_pairs++;
        #ifndef TRACK_CUTOFF_LOCATION
        std::cout << "\r" << "Progress: " << handled_pairs << "/" << n_pairs << std::flush;
        #endif
    }
    */

    #ifndef TRACK_CUTOFF_LOCATION
    std::cout << "\nDone!" << std::endl;
    #endif
    /*
    std::vector<convergent_pair> next = {};
    int iteration = 0;
    while (pairs.size() > 0) {
        next.reserve(pairs.size() * 2);
        for (auto pair: pairs) {
            if (!pair_meets_littlewood_criteria(pair)) {
                subdivide_and_append_new_pairs(pair, next);
            }
        }
        swap(pairs, next);
        next.clear();
        #ifndef TRACK_CUTOFF_LOCATION
        std::cout << "Pairs after iteration " << ++iteration << ": " << pairs.size() << std::endl;
        #endif
    }
    */
    
    #ifdef TRACK_MAX_DEPTH
    std::cout << "Max depth: " << max_depth << std::endl;
    #endif
    #ifdef TRACK_LARGEST
    std::cout << "Chonkiest: " << chonker << std::endl;
    #endif
}