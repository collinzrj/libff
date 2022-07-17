/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef DALEK_PAIRING_HPP_
#define DALEK_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/dalek/dalek_init.hpp>

namespace libff {

/* final exponentiation */

dalek_Fq6 dalek_final_exponentiation_last_chunk(const dalek_Fq6 &elt,
                                                    const dalek_Fq6 &elt_inv);
dalek_Fq6 dalek_final_exponentiation_first_chunk(const dalek_Fq6 &elt,
                                                     const dalek_Fq6 &elt_inv);
dalek_GT dalek_final_exponentiation(const dalek_Fq6 &elt);

/* Tate pairing */

struct dalek_Fq_conic_coefficients {
    dalek_Fq c_ZZ;
    dalek_Fq c_XY;
    dalek_Fq c_XZ;

    bool operator==(const dalek_Fq_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const dalek_Fq_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, dalek_Fq_conic_coefficients &cc);
};
typedef std::vector<dalek_Fq_conic_coefficients> dalek_tate_G1_precomp;

std::ostream& operator<<(std::ostream& out, const dalek_tate_G1_precomp &prec_P);
std::istream& operator>>(std::istream& in, dalek_tate_G1_precomp &prec_P);

struct dalek_tate_G2_precomp {
    dalek_Fq3 y0, eta;

    bool operator==(const dalek_tate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const dalek_tate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, dalek_tate_G2_precomp &prec_Q);
};

dalek_tate_G1_precomp dalek_tate_precompute_G1(const dalek_G1& P);
dalek_tate_G2_precomp dalek_tate_precompute_G2(const dalek_G2& Q);

dalek_Fq6 dalek_tate_miller_loop(const dalek_tate_G1_precomp &prec_P,
                                     const dalek_tate_G2_precomp &prec_Q);

dalek_Fq6 dalek_tate_pairing(const dalek_G1& P,
                                 const dalek_G2 &Q);
dalek_GT dalek_tate_reduced_pairing(const dalek_G1 &P,
                                        const dalek_G2 &Q);

/* ate pairing */

struct dalek_Fq3_conic_coefficients {
    dalek_Fq3 c_ZZ;
    dalek_Fq3 c_XY;
    dalek_Fq3 c_XZ;

    bool operator==(const dalek_Fq3_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const dalek_Fq3_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, dalek_Fq3_conic_coefficients &cc);
};
typedef std::vector<dalek_Fq3_conic_coefficients> dalek_ate_G2_precomp;

std::ostream& operator<<(std::ostream& out, const dalek_ate_G2_precomp &prec_Q);
std::istream& operator>>(std::istream& in, dalek_ate_G2_precomp &prec_Q);

struct dalek_ate_G1_precomp {
    dalek_Fq P_XY;
    dalek_Fq P_XZ;
    dalek_Fq P_ZZplusYZ;

    bool operator==(const dalek_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const dalek_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, dalek_ate_G1_precomp &prec_P);
};

dalek_ate_G1_precomp dalek_ate_precompute_G1(const dalek_G1& P);
dalek_ate_G2_precomp dalek_ate_precompute_G2(const dalek_G2& Q);

dalek_Fq6 dalek_ate_miller_loop(const dalek_ate_G1_precomp &prec_P,
                                    const dalek_ate_G2_precomp &prec_Q);
dalek_Fq6 dalek_ate_double_miller_loop(const dalek_ate_G1_precomp &prec_P1,
                                           const dalek_ate_G2_precomp &prec_Q1,
                                           const dalek_ate_G1_precomp &prec_P2,
                                           const dalek_ate_G2_precomp &prec_Q2);

dalek_Fq6 dalek_ate_pairing(const dalek_G1& P,
                                const dalek_G2 &Q);
dalek_GT dalek_ate_reduced_pairing(const dalek_G1 &P,
                                       const dalek_G2 &Q);

/* choice of pairing */

typedef dalek_ate_G1_precomp dalek_G1_precomp;
typedef dalek_ate_G2_precomp dalek_G2_precomp;

dalek_G1_precomp dalek_precompute_G1(const dalek_G1& P);
dalek_G2_precomp dalek_precompute_G2(const dalek_G2& Q);

dalek_Fq6 dalek_miller_loop(const dalek_G1_precomp &prec_P,
                                const dalek_G2_precomp &prec_Q);

dalek_Fq6 dalek_double_miller_loop(const dalek_G1_precomp &prec_P1,
                                       const dalek_G2_precomp &prec_Q1,
                                       const dalek_G1_precomp &prec_P2,
                                       const dalek_G2_precomp &prec_Q2);

dalek_Fq6 dalek_pairing(const dalek_G1& P,
                            const dalek_G2 &Q);

dalek_GT dalek_reduced_pairing(const dalek_G1 &P,
                                   const dalek_G2 &Q);

} // libff
#endif // DALEK_PAIRING_HPP_
