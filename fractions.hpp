/**
 * Copyright 2023 Topi Törmä, Matti Vapa
 */

#ifndef FRACTIONS
#define FRACTIONS

#include <vector>

namespace fractions {

    template <typename T>
    struct rational {
        T num;
        T den;
    };

    template <typename T>
    struct convergent {
        rational<T> current;
        rational<T> previous;
    };

    template <typename T>
    struct convergent_pair {
        convergent<T> alpha;
        convergent<T> beta;
    };

    template <typename T>
    convergent<T> next_convergent(const convergent<T>& base, T digit) {
        return {
            {
                digit * base.current.num + base.previous.num,
                digit * base.current.den + base.previous.den
            },
            base.current
        };
    }

    template <typename T>
    void subdivide(const convergent_pair<T>& pair, int N, std::vector<convergent_pair<T>>& results) {
        for (int i = 1; i < N; i++) {
            convergent next = next_convergent(pair.alpha, static_cast<T>(i));

            if (next.current.den < pair.beta.current.den) {
                results.push_back({next, pair.beta});
            } else {
                results.push_back({pair.beta, next});
            }
        }
    }

    template <typename T, typename F>
    void subdivide_with_cutoff_condition(const convergent_pair<T>& pair, int N, F&& cutoff_condition, std::vector<convergent_pair<T>>& results) {
        for (int i = 1; i < N; i++) {

            if (i > 1 && cutoff_condition(pair.alpha, pair.beta, i)) {
                break;
            }

            convergent next = next_convergent(pair.alpha, static_cast<T>(i));

            if (next.current.den < pair.beta.current.den) {
                results.push_back({next, pair.beta});
            } else {
                results.push_back({pair.beta, next});
            }
        }
    }

    template <typename T>
    std::vector<convergent_pair<T>> convergent_pairs(int N, uint subdivisions = 0) {
        
        std::vector<convergent<T>> candidates = {};
        // Create an initial set of candidate convergents.
        for (int i = 2; i < N; i++) {
            candidates.push_back({
                {static_cast<T>(1), static_cast<T>(i)},
                {static_cast<T>(0), static_cast<T>(1)}
            });
        }

        std::vector<convergent<T>> convergents = {};

        // "Refine" the initial set of candidates until we have a vector of convergents where
        // all of them have a sum of the denominators of the current and previous
        // iteration >= 2 * N.
        while (candidates.size() > 0) {
            // Pop a candidate from the end of the vector.
            auto candidate = std::move(candidates.back());
            candidates.pop_back();
            // If the sum of the denominators of the current and previous step of the of the candidate convergent
            // are equal or larger than 2 * N, add it to the vector of convergents.
            if (candidate.current.den + candidate.previous.den >= 2 * N) {
                convergents.push_back(candidate);
            } else {
                // Otherwise, we need to replace the candidate with the next iterations of the continued fraction
                // for all the following digits i where 1 <= i < N.
                for (int i = 1; i < N; i++) {
                    candidates.push_back(next_convergent(candidate, static_cast<T>(i)));
                }
            }
        }
        
        std::vector<convergent_pair<T>> pairs = {};
        
        // Make a list of all unique pairs of convergents.
        for (std::size_t i = 1; i < convergents.size(); i++) {
            for (std::size_t j = 0; j < i; j++) {
                convergent<T> a = convergents[i];
                convergent<T> b = convergents[j];
                

                // We can skip all the pairs where the condition (a_den/a_num) * (b_den/b_num) >= 2 * N,
                // as they will trivially fullfil the Littlewood criteria.
                T alpha_ratio = a.current.den / a.current.num;
                T beta_ratio = b.current.den / b.current.num;

                if (alpha_ratio * beta_ratio >= 2 * N) {
                    continue;
                }

                // Order the pair always so that the first convergent is the one with the smaller denominator.
                if (a.current.den < b.current.den) {
                    pairs.push_back({a, b});
                } else {
                    pairs.push_back({b, a});
                }
            }
        }

        // Optionally we can "subdivide" the pairs further (by replacing the pair with a list of pairs resulting
        // from replacing the first convergent in the pair with its following continued fraction iterations).
        // This has no mathematical bearing, but it allows us to split the computation over more CPU cores as
        // we'll begin the calcuation with a larger initial set of convergent pairs.
        if (subdivisions > 0) {
            for (uint i = 0; i < subdivisions; i++) {
                std::vector<convergent_pair<T>> next = {};
                for (auto& pair: pairs) {
                    subdivide(pair, N, next);
                }
                swap(next, pairs);
            }
        }
        
        return pairs;
    }
}

#endif