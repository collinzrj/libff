/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef EDWARDS_INIT_HPP_
#define EDWARDS_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp3.hpp>
#include <libff/algebra/fields/fp6_2over3.hpp>

namespace libff {

const mp_size_t dalek_r_bitcount = 254;
const mp_size_t dalek_q_bitcount = 254;

const mp_size_t dalek_r_limbs = (dalek_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t dalek_q_limbs = (dalek_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<dalek_r_limbs> dalek_modulus_r;
extern bigint<dalek_q_limbs> dalek_modulus_q;

typedef Fp_model<dalek_r_limbs, dalek_modulus_r> dalek_Fr;
typedef Fp_model<dalek_q_limbs, dalek_modulus_q> dalek_Fq;
typedef Fp3_model<dalek_q_limbs, dalek_modulus_q> dalek_Fq3;
typedef Fp6_2over3_model<dalek_q_limbs, dalek_modulus_q> dalek_Fq6;
typedef dalek_Fq6 dalek_GT;

// parameters for Edwards curve E_{1,d}(F_q)
extern dalek_Fq dalek_coeff_a;
extern dalek_Fq dalek_coeff_d;
// parameters for twisted Edwards curve E_{a',d'}(F_q^3)
extern dalek_Fq3 dalek_twist;
extern dalek_Fq3 dalek_twist_coeff_a;
extern dalek_Fq3 dalek_twist_coeff_d;
extern dalek_Fq dalek_twist_mul_by_a_c0;
extern dalek_Fq dalek_twist_mul_by_a_c1;
extern dalek_Fq dalek_twist_mul_by_a_c2;
extern dalek_Fq dalek_twist_mul_by_d_c0;
extern dalek_Fq dalek_twist_mul_by_d_c1;
extern dalek_Fq dalek_twist_mul_by_d_c2;
extern dalek_Fq dalek_twist_mul_by_q_Y;
extern dalek_Fq dalek_twist_mul_by_q_Z;

// parameters for pairing
extern bigint<dalek_q_limbs> dalek_ate_loop_count;
extern bigint<6*dalek_q_limbs> dalek_final_exponent;
extern bigint<dalek_q_limbs> dalek_final_exponent_last_chunk_abs_of_w0;
extern bool dalek_final_exponent_last_chunk_is_w0_neg;
extern bigint<dalek_q_limbs> dalek_final_exponent_last_chunk_w1;

void init_dalek_params();

class dalek_G1;
class dalek_G2;

} // libff
#endif // EDWARDS_INIT_HPP_
