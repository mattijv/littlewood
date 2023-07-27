/**
 * Copyright 2023 Topi Törmä, Matti Vapa
 */

#ifndef LITTLEWOOD
#define LITTLEWOOD

#include <algorithm>
#include "fractions.hpp"
#include "modular_math.hpp"

namespace LW {
    template <typename T>
    T littlewood(T Q, T a, T b, int N) {
        T result = 2 * N * Q * (Q + a) * (Q + b);
        return result;
    }

    template <typename T>
    struct littlewood_result {
        T best_q;
        bool meets_criteria;
    };
    
    template <typename T>
    littlewood_result<T> meets_littlewood_criteria(const fractions::convergent_pair<T>& pair, int N) {
        
        fractions::convergent<T> alpha = pair.alpha;
        fractions::convergent<T> beta = pair.beta;

        T alpha_sum = alpha.current.den + alpha.previous.den;
        T beta_sum = beta.current.den + beta.previous.den;

        T epsilon = alpha.current.den * alpha_sum * beta.current.den * beta_sum;
        
        // Quickly check if the bigger convergent denominator is good.
        // Corresponds to step 2 (a) of the algorithm in the paper.
        T alpha_remainder = modular_math::remainder_with_least_absolute_value(beta.current.den, alpha.current);
        T littlewood_quantity = littlewood(beta.current.den, alpha_remainder * alpha_sum, static_cast<T>(0), N);
        if (littlewood_quantity < epsilon) {
            return {static_cast<T>(0), true};
        }

        // Track the Q that gives the lowest value for the "Littlewood quantity".
        T best_q = beta.current.den;
        T lowest_littlewood_quantity = littlewood_quantity;

        T alpha_diff = alpha.current.den - alpha.previous.den;
        T beta_diff = beta.current.den - beta.previous.den;

        T denominators[4] = {
            alpha.previous.den,
            alpha_diff,
            beta.previous.den,
            beta_diff
        };

        T alpha_den_beta_num = (alpha.previous.den * beta.current.num) % beta.current.den;
        T alpha_den_comp_beta_num = (alpha_diff * beta.current.num) % beta.current.den;
        T beta_den_alpha_num = (beta.previous.den * alpha.current.num) % alpha.current.den;
        T beta_den_comp_alpha_num = (beta_diff * alpha.current.num) % alpha.current.den;
        T multiple_alpha = (alpha.current.den * beta.current.num) % beta.current.den;

        T ab_rem = alpha_den_beta_num;
        T acb_rem = alpha_den_comp_beta_num;
        T ba_rem = beta_den_alpha_num;
        T bca_rem = beta_den_comp_alpha_num;
        T ma_rem = multiple_alpha;

        // max_remainder corresponds to the "suitably large integer" described in step 2 (b) of the algorithm in the article.
        // Specifically this corresponds to M + 1 to make the condition of the while loop below a bit neater.
        T max_remainder = static_cast<T>(1) + std::max(static_cast<T>(N), pair.beta.previous.den / (2 * N * N));
        // target_remainder corresponds to the r described in step 2 (b) of the algorithm in the article.
        T target_remainder = static_cast<T>(1);
        
        // Check all r, where 1 <= r <= M
        while (target_remainder < max_remainder) {

            // The following checks correspond to the substeps i - v in the step 2 (b) of the algorithm.
            littlewood_quantity = littlewood(denominators[0], target_remainder * alpha_sum, std::min(ab_rem, beta.current.den - ab_rem) * beta_sum, N);
            if (littlewood_quantity < epsilon) {
                return {static_cast<T>(0), true};
            } else if (littlewood_quantity < lowest_littlewood_quantity) {
                lowest_littlewood_quantity = littlewood_quantity;
                best_q = denominators[0];
            }

            littlewood_quantity = littlewood(denominators[1], target_remainder * alpha_sum, std::min(acb_rem, beta.current.den - acb_rem) * beta_sum, N);
            if (littlewood_quantity < epsilon) {
                return {static_cast<T>(0), true};
            } else if (littlewood_quantity < lowest_littlewood_quantity) {
                lowest_littlewood_quantity = littlewood_quantity;
                best_q = denominators[1];
            }

            littlewood_quantity = littlewood(denominators[2], std::min(ba_rem, alpha.current.den - ba_rem) * alpha_sum, target_remainder * beta_sum, N);
            if (littlewood_quantity < epsilon) {
                return {static_cast<T>(0), true};
            } else if (littlewood_quantity < lowest_littlewood_quantity) {
                lowest_littlewood_quantity = littlewood_quantity;
                best_q = denominators[2];
            }

            littlewood_quantity = littlewood(denominators[3], std::min(bca_rem, alpha.current.den - bca_rem) * alpha_sum, target_remainder * beta_sum, N);
            if (littlewood_quantity < epsilon) {
                return {static_cast<T>(0), true};
            } else if (littlewood_quantity < lowest_littlewood_quantity) {
                lowest_littlewood_quantity = littlewood_quantity;
                best_q = denominators[3];
            }

            littlewood_quantity = littlewood(target_remainder * alpha.current.den, static_cast<T>(0), std::min(ma_rem, beta.current.den - ma_rem) * beta_sum, N);
            if (littlewood_quantity < epsilon) {
                return {static_cast<T>(0), true};
            } else if (littlewood_quantity < lowest_littlewood_quantity) {
                lowest_littlewood_quantity = littlewood_quantity;
                best_q = target_remainder * alpha.current.den;
            }

            target_remainder++;

            // The following section uses recursive properties to minimize the need to perform modulo operations.

            T ab_sum = denominators[0] + alpha.previous.den;

            if (ab_sum <= alpha.current.den) {
                denominators[0] = ab_sum;
                ab_rem = modular_math::modular_addition(ab_rem, alpha_den_beta_num, beta.current.den);
            } else {
                denominators[0] = ab_sum - alpha.current.den;
                ab_rem = modular_math::modular_substraction(ab_rem, alpha_den_comp_beta_num, beta.current.den);
            }

            T acb_sum = denominators[1] + alpha_diff;
            if (acb_sum <= alpha.current.den) {
                denominators[1] = acb_sum;
                acb_rem = modular_math::modular_addition(acb_rem, alpha_den_comp_beta_num, beta.current.den);
            }
            else {
                denominators[1] = acb_sum - alpha.current.den;
                acb_rem = modular_math::modular_substraction(acb_rem, alpha_den_beta_num, beta.current.den);
            }

            T ba_sum = denominators[2] + beta.previous.den;
            if (ba_sum <= beta.current.den) {
                denominators[2] = ba_sum;
                ba_rem = modular_math::modular_addition(ba_rem, beta_den_alpha_num, alpha.current.den);
            }
            else {
                denominators[2] = ba_sum - beta.current.den;
                ba_rem = modular_math::modular_substraction(ba_rem, beta_den_comp_alpha_num, alpha.current.den);
            }

            T bca_sum = denominators[3] + beta_diff;
            if (bca_sum <= beta.current.den) {
                denominators[3] = bca_sum;
                bca_rem = modular_math::modular_addition(bca_rem, beta_den_comp_alpha_num, alpha.current.den);
            }
            else {
                denominators[3] = bca_sum - beta.current.den;
                bca_rem = modular_math::modular_substraction(bca_rem, beta_den_alpha_num, alpha.current.den);
            }

            ma_rem = modular_math::modular_addition(ma_rem, multiple_alpha, beta.current.den);
        }
        
        // The current pair did not meet the Littlewood criteria for any value of r, so we just return the best q.
        return {best_q, false};
    }

    /**
     * The cutoff method implements the check using equation (31) of the article
     */
    template <typename T>
    bool littlewood_cutoff_reached(T best_q, const fractions::convergent<T> alpha, const fractions::convergent<T> beta, int next_digit, int N) {
        T new_alpha_sum = next_digit * alpha.current.den + alpha.previous.den;
        T beta_sum = beta.current.den + beta.previous.den;
        T new_epsilon = alpha.current.den * new_alpha_sum * beta.current.den * beta_sum;
        T a = new_alpha_sum * modular_math::remainder_with_least_absolute_value(best_q, alpha.current);
        T b = beta_sum * modular_math::remainder_with_least_absolute_value(best_q, beta.current);
        return littlewood(best_q, a, b, N) < new_epsilon;
    }
}

#endif