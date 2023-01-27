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
        convergent<T> primary;
        convergent<T> secondary;

        if (pair.alpha.current.den < pair.beta.current.den) {
            primary = pair.alpha;
            secondary = pair.beta;
        } else {
            primary = pair.beta;
            secondary = pair.alpha;
        }

        for (int i = 1; i < N; i++) {
            convergent next = next_convergent(primary, static_cast<T>(i));
            
            if (next.current.den < secondary.current.den) {
                results.push_back({next, secondary});
            } else {
                results.push_back({secondary, next});
            }
        }
    }

    template <typename T>
    std::vector<convergent_pair<T>> convergent_pairs(int N, uint subdivisions = 0) {
        std::vector<convergent<T>> convergents = {};
        
        for (std::size_t i = 2; i < N; i++) {
            convergents.push_back({
                {static_cast<T>(1), static_cast<T>(i)},
                {static_cast<T>(0), static_cast<T>(1)}
            });
        }
        
        std::vector<convergent_pair<T>> pairs = {};
        
        for (std::size_t i = 0; i < convergents.size(); i++) {
            for (std::size_t j = 0; j <= i; j++) {
                convergent<T> a = convergents[i];
                convergent<T> b = convergents[j];
                
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

        if (subdivisions > 0) {
            for (int i = 0; i < subdivisions; i++) {
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