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
    bool meets_littlewood_criteria(const fractions::convergent_pair<T>& pair, int N) {
        
        fractions::convergent<T> alpha = pair.alpha;
        fractions::convergent<T> beta = pair.beta;

        T alpha_sum = alpha.current.den + alpha.previous.den;
        T beta_sum = beta.current.den + beta.previous.den;

        T epsilon = alpha.current.den * alpha_sum * beta.current.den * beta_sum;
        
        //Quickly check if the bigger convergent denominator is good
        T alpha_remainder = modular_math::remainder_with_least_absolute_value(beta.current.den, alpha.current);
        if (littlewood(beta.current.den, static_cast<T>(0), alpha_remainder * (alpha.current.den + alpha.previous.den), N) < epsilon) {
            return true;
        }

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

        T max_remainder = 1 + std::max(static_cast<T>(N), pair.beta.previous.den / (N * N));
        T target_remainder = static_cast<T>(1);
        
        while (target_remainder < max_remainder) {

            if (littlewood(denominators[0], target_remainder * alpha_sum, std::min(ab_rem, beta.current.den - ab_rem) * beta_sum, N) < epsilon) {
                return true;
            }

            if (littlewood(denominators[1], target_remainder * alpha_sum, std::min(acb_rem, beta.current.den - acb_rem) * beta_sum, N) < epsilon) {
                return true;
            }

            if (littlewood(denominators[2], std::min(ba_rem, alpha.current.den - ba_rem) * alpha_sum, target_remainder * beta_sum, N) < epsilon) {
                return true;
            }

            if (littlewood(denominators[3], std::min(bca_rem, alpha.current.den - bca_rem) * alpha_sum, target_remainder * beta_sum, N) < epsilon) {
                return true;
            }

            if (littlewood(target_remainder * alpha.current.den, static_cast<T>(0), std::min(ma_rem, beta.current.den - ma_rem) * beta_sum, N) < epsilon) {
                return true;
            }

            target_remainder++;

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
        
        return false;
    }
    
}

#endif