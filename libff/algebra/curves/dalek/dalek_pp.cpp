/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/dalek/dalek_pp.hpp>

namespace libff {

void dalek_pp::init_public_params()
{
    init_dalek_params();
}

dalek_GT dalek_pp::final_exponentiation(const dalek_Fq6 &elt)
{
    return dalek_final_exponentiation(elt);
}

dalek_G1_precomp dalek_pp::precompute_G1(const dalek_G1 &P)
{
    return dalek_precompute_G1(P);
}

dalek_G2_precomp dalek_pp::precompute_G2(const dalek_G2 &Q)
{
    return dalek_precompute_G2(Q);
}

dalek_Fq6 dalek_pp::miller_loop(const dalek_G1_precomp &prec_P,
                                    const dalek_G2_precomp &prec_Q)
{
    return dalek_miller_loop(prec_P, prec_Q);
}

dalek_Fq6 dalek_pp::double_miller_loop(const dalek_G1_precomp &prec_P1,
                                           const dalek_G2_precomp &prec_Q1,
                                           const dalek_G1_precomp &prec_P2,
                                           const dalek_G2_precomp &prec_Q2)
{
    return dalek_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

dalek_Fq6 dalek_pp::pairing(const dalek_G1 &P,
                                const dalek_G2 &Q)
{
    return dalek_pairing(P, Q);
}

dalek_Fq6 dalek_pp::reduced_pairing(const dalek_G1 &P,
                                        const dalek_G2 &Q)
{
    return dalek_reduced_pairing(P, Q);
}

} // libff
