/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef DALEK_PP_HPP_
#define DALEK_PP_HPP_
#include <libff/algebra/curves/dalek/dalek_g1.hpp>
#include <libff/algebra/curves/dalek/dalek_g2.hpp>
#include <libff/algebra/curves/dalek/dalek_init.hpp>
#include <libff/algebra/curves/dalek/dalek_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class dalek_pp {
public:
    typedef dalek_Fr Fp_type;
    typedef dalek_G1 G1_type;
    typedef dalek_G2 G2_type;
    typedef dalek_G1_precomp G1_precomp_type;
    typedef dalek_G2_precomp G2_precomp_type;
    typedef dalek_Fq Fq_type;
    typedef dalek_Fq3 Fqe_type;
    typedef dalek_Fq6 Fqk_type;
    typedef dalek_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static dalek_GT final_exponentiation(const dalek_Fq6 &elt);
    static dalek_G1_precomp precompute_G1(const dalek_G1 &P);
    static dalek_G2_precomp precompute_G2(const dalek_G2 &Q);
    static dalek_Fq6 miller_loop(const dalek_G1_precomp &prec_P,
                                   const dalek_G2_precomp &prec_Q);
    static dalek_Fq6 double_miller_loop(const dalek_G1_precomp &prec_P1,
                                          const dalek_G2_precomp &prec_Q1,
                                          const dalek_G1_precomp &prec_P2,
                                          const dalek_G2_precomp &prec_Q2);
    /* the following are used in test files */
    static dalek_Fq6 pairing(const dalek_G1 &P,
                               const dalek_G2 &Q);
    static dalek_Fq6 reduced_pairing(const dalek_G1 &P,
                                       const dalek_G2 &Q);
};

} // libff
#endif // DALEK_PP_HPP_
