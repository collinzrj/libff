/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <cassert>

#include <libff/algebra/curves/dalek/dalek_g1.hpp>
#include <libff/algebra/curves/dalek/dalek_g2.hpp>
#include <libff/algebra/curves/dalek/dalek_init.hpp>
#include <libff/algebra/curves/dalek/dalek_pairing.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

bool dalek_Fq_conic_coefficients::operator==(const dalek_Fq_conic_coefficients &other) const
{
    return (this->c_ZZ == other.c_ZZ &&
            this->c_XY == other.c_XY &&
            this->c_XZ == other.c_XZ);
}

std::ostream& operator<<(std::ostream &out, const dalek_Fq_conic_coefficients &cc)
{
    out << cc.c_ZZ << OUTPUT_SEPARATOR << cc.c_XY << OUTPUT_SEPARATOR << cc.c_XZ;
    return out;
}

std::istream& operator>>(std::istream &in, dalek_Fq_conic_coefficients &cc)
{
    in >> cc.c_ZZ;
    consume_OUTPUT_SEPARATOR(in);
    in >> cc.c_XY;
    consume_OUTPUT_SEPARATOR(in);
    in >> cc.c_XZ;
    return in;
}

std::ostream& operator<<(std::ostream& out, const dalek_tate_G1_precomp &prec_P)
{
    out << prec_P.size() << "\n";
    for (const dalek_Fq_conic_coefficients &cc : prec_P)
    {
        out << cc << OUTPUT_NEWLINE;
    }

    return out;
}

std::istream& operator>>(std::istream& in, dalek_tate_G1_precomp &prec_P)
{
    prec_P.clear();

    size_t s;
    in >> s;

    consume_newline(in);
    prec_P.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        dalek_Fq_conic_coefficients cc;
        in >> cc;
        consume_OUTPUT_NEWLINE(in);
        prec_P.emplace_back(cc);
    }

    return in;
}

bool dalek_tate_G2_precomp::operator==(const dalek_tate_G2_precomp &other) const
{
    return (this->y0 == other.y0 &&
            this->eta == other.eta);
}

std::ostream& operator<<(std::ostream &out, const dalek_tate_G2_precomp &prec_Q)
{
    out << prec_Q.y0 << OUTPUT_SEPARATOR << prec_Q.eta;
    return out;
}

std::istream& operator>>(std::istream &in, dalek_tate_G2_precomp &prec_Q)
{
    in >> prec_Q.y0;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.eta;
    return in;
}

bool dalek_Fq3_conic_coefficients::operator==(const dalek_Fq3_conic_coefficients &other) const
{
    return (this->c_ZZ == other.c_ZZ &&
            this->c_XY == other.c_XY &&
            this->c_XZ == other.c_XZ);
}

std::ostream& operator<<(std::ostream &out, const dalek_Fq3_conic_coefficients &cc)
{
    out << cc.c_ZZ << OUTPUT_SEPARATOR << cc.c_XY << OUTPUT_SEPARATOR << cc.c_XZ;
    return out;
}

std::istream& operator>>(std::istream &in, dalek_Fq3_conic_coefficients &cc)
{
    in >> cc.c_ZZ;
    consume_OUTPUT_SEPARATOR(in);
    in >> cc.c_XY;
    consume_OUTPUT_SEPARATOR(in);
    in >> cc.c_XZ;
    return in;
}

std::ostream& operator<<(std::ostream& out, const dalek_ate_G2_precomp &prec_Q)
{
    out << prec_Q.size() << "\n";
    for (const dalek_Fq3_conic_coefficients &cc : prec_Q)
    {
        out << cc << OUTPUT_NEWLINE;
    }

    return out;
}

std::istream& operator>>(std::istream& in, dalek_ate_G2_precomp &prec_Q)
{
    prec_Q.clear();

    size_t s;
    in >> s;

    consume_newline(in);

    prec_Q.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        dalek_Fq3_conic_coefficients cc;
        in >> cc;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.emplace_back(cc);
    }

    return in;
}

bool dalek_ate_G1_precomp::operator==(const dalek_ate_G1_precomp &other) const
{
    return (this->P_XY == other.P_XY &&
            this->P_XZ == other.P_XZ &&
            this->P_ZZplusYZ == other.P_ZZplusYZ);
}

std::ostream& operator<<(std::ostream &out, const dalek_ate_G1_precomp &prec_P)
{
    out << prec_P.P_XY << OUTPUT_SEPARATOR << prec_P.P_XZ << OUTPUT_SEPARATOR << prec_P.P_ZZplusYZ;

    return out;
}

std::istream& operator>>(std::istream &in, dalek_ate_G1_precomp &prec_P)
{
    in >> prec_P.P_XY >> prec_P.P_XZ >> prec_P.P_ZZplusYZ;

    return in;
}

/* final exponentiations */
dalek_Fq6 dalek_final_exponentiation_last_chunk(const dalek_Fq6 &elt, const dalek_Fq6 &elt_inv)
{
    enter_block("Call to dalek_final_exponentiation_last_chunk");
    const dalek_Fq6 elt_q = elt.Frobenius_map(1);
    dalek_Fq6 w1_part = elt_q.cyclotomic_exp(dalek_final_exponent_last_chunk_w1);
    dalek_Fq6 w0_part;
    if (dalek_final_exponent_last_chunk_is_w0_neg)
    {
    	w0_part = elt_inv.cyclotomic_exp(dalek_final_exponent_last_chunk_abs_of_w0);
    } else {
    	w0_part = elt.cyclotomic_exp(dalek_final_exponent_last_chunk_abs_of_w0);
    }
    dalek_Fq6 result = w1_part * w0_part;
    leave_block("Call to dalek_final_exponentiation_last_chunk");

    return result;
}

dalek_Fq6 dalek_final_exponentiation_first_chunk(const dalek_Fq6 &elt, const dalek_Fq6 &elt_inv)
{
    enter_block("Call to dalek_final_exponentiation_first_chunk");

    /* (q^3-1)*(q+1) */

    /* elt_q3 = elt^(q^3) */
    const dalek_Fq6 elt_q3 = elt.Frobenius_map(3);
    /* elt_q3_over_elt = elt^(q^3-1) */
    const dalek_Fq6 elt_q3_over_elt = elt_q3 * elt_inv;
    /* alpha = elt^((q^3-1) * q) */
    const dalek_Fq6 alpha = elt_q3_over_elt.Frobenius_map(1);
    /* beta = elt^((q^3-1)*(q+1) */
    const dalek_Fq6 beta = alpha * elt_q3_over_elt;
    leave_block("Call to dalek_final_exponentiation_first_chunk");
    return beta;
}

dalek_GT dalek_final_exponentiation(const dalek_Fq6 &elt)
{
    enter_block("Call to dalek_final_exponentiation");
    const dalek_Fq6 elt_inv = elt.inverse();
    const dalek_Fq6 elt_to_first_chunk = dalek_final_exponentiation_first_chunk(elt, elt_inv);
    const dalek_Fq6 elt_inv_to_first_chunk = dalek_final_exponentiation_first_chunk(elt_inv, elt);
    dalek_GT result = dalek_final_exponentiation_last_chunk(elt_to_first_chunk, elt_inv_to_first_chunk);
    leave_block("Call to dalek_final_exponentiation");

    return result;
}

dalek_tate_G2_precomp dalek_tate_precompute_G2(const dalek_G2& Q)
{
    enter_block("Call to dalek_tate_precompute_G2");
    dalek_G2 Qcopy = Q;
    Qcopy.to_affine_coordinates();
    dalek_tate_G2_precomp result;
    result.y0 = Qcopy.Y * Qcopy.Z.inverse(); // Y/Z
    result.eta = (Qcopy.Z+Qcopy.Y) * dalek_Fq6::mul_by_non_residue(Qcopy.X).inverse(); // (Z+Y)/(nqr*X)
    leave_block("Call to dalek_tate_precompute_G2");

    return result;
}

struct extended_dalek_G1_projective {
    dalek_Fq X;
    dalek_Fq Y;
    dalek_Fq Z;
    dalek_Fq T;

    void print() const
        {
            printf("extended dalek_G1 projective X/Y/Z/T:\n");
            X.print();
            Y.print();
            Z.print();
            T.print();
        }

    void test_invariant() const
        {
            assert(T*Z == X*Y);
        }
};

void doubling_step_for_miller_loop(extended_dalek_G1_projective &current,
                                   dalek_Fq_conic_coefficients &cc)
{
    const dalek_Fq &X = current.X, &Y = current.Y, &Z = current.Z, &T = current.T;
    const dalek_Fq A = X.squared();     // A    = X1^2
    const dalek_Fq B = Y.squared();     // B    = Y1^2
    const dalek_Fq C = Z.squared();     // C    = Z1^2
    const dalek_Fq D = (X+Y).squared(); // D    = (X1+Y1)^2
    const dalek_Fq E = (Y+Z).squared(); // E    = (Y1+Z1)^2
    const dalek_Fq F = D-(A+B);         // F    = D-(A+B)
    const dalek_Fq G = E-(B+C);         // G    = E-(B+C)
    const dalek_Fq &H = A;              // H    = A (dalek_a=1)
    const dalek_Fq I = H+B;             // I    = H+B
    const dalek_Fq J = C-I;             // J    = C-I
    const dalek_Fq K = J+C;             // K    = J+C

    cc.c_ZZ = Y*(T-X);            // c_ZZ = 2*Y1*(T1-X1)
    cc.c_ZZ = cc.c_ZZ + cc.c_ZZ;

    cc.c_XY = J+J+G;              // c_XY = 2*J+G
    cc.c_XZ = X*T-B;              // c_XZ = 2*(X1*T1-B) (dalek_a=1)
    cc.c_XZ = cc.c_XZ + cc.c_XZ;

    current.X = F*K;              // X3 = F*K
    current.Y = I*(B-H);          // Y3 = I*(B-H)
    current.Z = I*K;              // Z3 = I*K
    current.T = F*(B-H);          // T3 = F*(B-H)

#ifdef DEBUG
    current.test_invariant();
#endif
}

void full_addition_step_for_miller_loop(const extended_dalek_G1_projective &base,
                                        extended_dalek_G1_projective &current,
                                        dalek_Fq_conic_coefficients &cc)
{
    const dalek_Fq &X1 = current.X, &Y1 = current.Y, &Z1 = current.Z, &T1 = current.T;
    const dalek_Fq &X2 = base.X, &Y2 =  base.Y, &Z2 = base.Z, &T2 = base.T;

    const dalek_Fq A = X1*X2;               // A    = X1*X2
    const dalek_Fq B = Y1*Y2;               // B    = Y1*Y2
    const dalek_Fq C = Z1*T2;               // C    = Z1*T2
    const dalek_Fq D = T1*Z2;               // D    = T1*Z2
    const dalek_Fq E = D+C;                 // E    = D+C
    const dalek_Fq F = (X1-Y1)*(X2+Y2)+B-A; // F    = (X1-Y1)*(X2+Y2)+B-A
    const dalek_Fq G = B + A;               // G    = B + A (dalek_a=1)
    const dalek_Fq H = D-C;                 // H    = D-C
    const dalek_Fq I = T1*T2;               // I    = T1*T2

    cc.c_ZZ = (T1-X1)*(T2+X2)-I+A;    // c_ZZ = (T1-X1)*(T2+X2)-I+A
    cc.c_XY = X1*Z2-X2*Z1+F;          // c_XY = X1*Z2-X2*Z1+F
    cc.c_XZ = (Y1-T1)*(Y2+T2)-B+I-H;  // c_XZ = (Y1-T1)*(Y2+T2)-B+I-H
    current.X = E*F;                  // X3   = E*F
    current.Y = G*H;                  // Y3   = G*H
    current.Z = F*G;                  // Z3   = F*G
    current.T = E*H;                  // T3   = E*H

#ifdef DEBUG
    current.test_invariant();
#endif
}

void mixed_addition_step_for_miller_loop(const extended_dalek_G1_projective &base,
                                         extended_dalek_G1_projective &current,
                                         dalek_Fq_conic_coefficients &cc)
{
    const dalek_Fq &X1 = current.X, &Y1 = current.Y, &Z1 = current.Z, &T1 = current.T;
    const dalek_Fq &X2 = base.X, &Y2 =  base.Y, &T2 = base.T;

    const dalek_Fq A = X1*X2;               // A    = X1*X2
    const dalek_Fq B = Y1*Y2;               // B    = Y1*Y2
    const dalek_Fq C = Z1*T2;               // C    = Z1*T2
    const dalek_Fq D = T1;                  // D    = T1*Z2
    const dalek_Fq E = D+C;                 // E    = D+C
    const dalek_Fq F = (X1-Y1)*(X2+Y2)+B-A; // F    = (X1-Y1)*(X2+Y2)+B-A
    const dalek_Fq G = B + A;               // G    = B + A (dalek_a=1)
    const dalek_Fq H = D-C;                 // H    = D-C
    const dalek_Fq I = T1*T2;               // I    = T1*T2

    cc.c_ZZ = (T1-X1)*(T2+X2)-I+A;    // c_ZZ = (T1-X1)*(T2+X2)-I+A
    cc.c_XY = X1-X2*Z1+F;             // c_XY = X1*Z2-X2*Z1+F
    cc.c_XZ = (Y1-T1)*(Y2+T2)-B+I-H;  // c_XZ = (Y1-T1)*(Y2+T2)-B+I-H
    current.X = E*F;                  // X3   = E*F
    current.Y = G*H;                  // Y3   = G*H
    current.Z = F*G;                  // Z3   = F*G
    current.T = E*H;                  // T3   = E*H

#ifdef DEBUG
    current.test_invariant();
#endif
}

dalek_tate_G1_precomp dalek_tate_precompute_G1(const dalek_G1& P)
{
    enter_block("Call to dalek_tate_precompute_G1");
    dalek_tate_G1_precomp result;

    dalek_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();

    extended_dalek_G1_projective P_ext;
    P_ext.X = Pcopy.X;
    P_ext.Y = Pcopy.Y;
    P_ext.Z = Pcopy.Z;
    P_ext.T = Pcopy.X*Pcopy.Y;

    extended_dalek_G1_projective R = P_ext;

    bool found_one = false;
    for (long i = dalek_modulus_r.max_bits(); i >= 0; --i)
    {
        const bool bit = dalek_modulus_r.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           dalek_modulus_r (skipping leading zeros) in MSB to LSB
           order */
        dalek_Fq_conic_coefficients cc;
        doubling_step_for_miller_loop(R, cc);
        result.push_back(cc);

        if (bit)
        {
            mixed_addition_step_for_miller_loop(P_ext, R, cc);
            result.push_back(cc);
        }
    }

    leave_block("Call to dalek_tate_precompute_G1");
    return result;
}

dalek_Fq6 dalek_tate_miller_loop(const dalek_tate_G1_precomp &prec_P,
                          const dalek_tate_G2_precomp &prec_Q)
{
    enter_block("Call to dalek_tate_miller_loop");

    dalek_Fq6 f = dalek_Fq6::one();

    bool found_one = false;
    size_t idx = 0;
    for (long i = dalek_modulus_r.max_bits()-1; i >= 0; --i)
    {
        const bool bit = dalek_modulus_r.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           dalek_modulus_r (skipping leading zeros) in MSB to LSB
           order */
        dalek_Fq_conic_coefficients cc = prec_P[idx++];
        dalek_Fq6 g_RR_at_Q = dalek_Fq6(dalek_Fq3(cc.c_XZ, dalek_Fq(0l), dalek_Fq(0l)) + cc.c_XY * prec_Q.y0,
                                            cc.c_ZZ * prec_Q.eta);
        f = f.squared() * g_RR_at_Q;
        if (bit)
        {
            cc = prec_P[idx++];

            dalek_Fq6 g_RP_at_Q = dalek_Fq6(dalek_Fq3(cc.c_XZ, dalek_Fq(0l), dalek_Fq(0l)) + cc.c_XY * prec_Q.y0,
                                                cc.c_ZZ * prec_Q.eta);
            f = f * g_RP_at_Q;
        }
    }
    leave_block("Call to dalek_tate_miller_loop");

    return f;
}

dalek_Fq6 dalek_tate_pairing(const dalek_G1& P, const dalek_G2 &Q)
{
    enter_block("Call to dalek_tate_pairing");
    dalek_tate_G1_precomp prec_P = dalek_tate_precompute_G1(P);
    dalek_tate_G2_precomp prec_Q = dalek_tate_precompute_G2(Q);
    dalek_Fq6 result = dalek_tate_miller_loop(prec_P, prec_Q);
    leave_block("Call to dalek_tate_pairing");
    return result;
}

dalek_GT dalek_tate_reduced_pairing(const dalek_G1 &P, const dalek_G2 &Q)
{
    enter_block("Call to dalek_tate_reduced_pairing");
    const dalek_Fq6 f = dalek_tate_pairing(P, Q);
    const dalek_GT result = dalek_final_exponentiation(f);
    leave_block("Call to dalek_tate_reduce_pairing");
    return result;
}

struct extended_dalek_G2_projective {
    dalek_Fq3 X;
    dalek_Fq3 Y;
    dalek_Fq3 Z;
    dalek_Fq3 T;

    void print() const
        {
            printf("extended dalek_G2 projective X/Y/Z/T:\n");
            X.print();
            Y.print();
            Z.print();
            T.print();
        }

    void test_invariant() const
        {
            assert(T*Z == X*Y);
        }
};

void doubling_step_for_flipped_miller_loop(extended_dalek_G2_projective &current,
                                           dalek_Fq3_conic_coefficients &cc)
{
    const dalek_Fq3 &X = current.X, &Y = current.Y, &Z = current.Z, &T = current.T;
    const dalek_Fq3 A = X.squared();     // A    = X1^2
    const dalek_Fq3 B = Y.squared();     // B    = Y1^2
    const dalek_Fq3 C = Z.squared();     // C    = Z1^2
    const dalek_Fq3 D = (X+Y).squared(); // D    = (X1+Y1)^2
    const dalek_Fq3 E = (Y+Z).squared(); // E    = (Y1+Z1)^2
    const dalek_Fq3 F = D-(A+B);         // F    = D-(A+B)
    const dalek_Fq3 G = E-(B+C);         // G    = E-(B+C)
    const dalek_Fq3 H = dalek_G2::mul_by_a(A); // dalek_param_twist_coeff_a is 1 * X for us
                                   // H    = twisted_a * A
    const dalek_Fq3 I = H+B;             // I    = H+B
    const dalek_Fq3 J = C-I;             // J    = C-I
    const dalek_Fq3 K = J+C;             // K    = J+C

    cc.c_ZZ = Y*(T-X);            // c_ZZ = 2*Y1*(T1-X1)
    cc.c_ZZ = cc.c_ZZ + cc.c_ZZ;

    // c_XY = 2*(C-dalek_a * A * delta_3-B)+G (dalek_a = 1 for us)
    cc.c_XY = C - dalek_G2::mul_by_a(A) - B; // dalek_param_twist_coeff_a is 1 * X for us
    cc.c_XY = cc.c_XY + cc.c_XY + G;

    // c_XZ = 2*(dalek_a*X1*T1*delta_3-B) (dalek_a = 1 for us)
    cc.c_XZ = dalek_G2::mul_by_a(X * T) - B; // dalek_param_twist_coeff_a is 1 * X for us
    cc.c_XZ = cc.c_XZ + cc.c_XZ;

    current.X = F*K;              // X3 = F*K
    current.Y = I*(B-H);          // Y3 = I*(B-H)
    current.Z = I*K;              // Z3 = I*K
    current.T = F*(B-H);          // T3 = F*(B-H)
#ifdef DEBUG
    current.test_invariant();
#endif
}

void full_addition_step_for_flipped_miller_loop(const extended_dalek_G2_projective &base,
                                                extended_dalek_G2_projective &current,
                                                dalek_Fq3_conic_coefficients &cc)
{
    const dalek_Fq3 &X1 = current.X, &Y1 = current.Y, &Z1 = current.Z, &T1 = current.T;
    const dalek_Fq3 &X2 = base.X, &Y2 =  base.Y, &Z2 = base.Z, &T2 = base.T;

    const dalek_Fq3 A = X1*X2;               // A    = X1*X2
    const dalek_Fq3 B = Y1*Y2;               // B    = Y1*Y2
    const dalek_Fq3 C = Z1*T2;               // C    = Z1*T2
    const dalek_Fq3 D = T1*Z2;               // D    = T1*Z2
    const dalek_Fq3 E = D+C;                 // E    = D+C
    const dalek_Fq3 F = (X1-Y1)*(X2+Y2)+B-A; // F    = (X1-Y1)*(X2+Y2)+B-A
    // G = B + twisted_dalek_a * A
    const dalek_Fq3 G = B + dalek_G2::mul_by_a(A); // dalek_param_twist_coeff_a is 1*X for us
    const dalek_Fq3 H = D-C;                 // H    = D-C
    const dalek_Fq3 I = T1*T2;               // I    = T1*T2

    // c_ZZ = delta_3* ((T1-X1)*(T2+X2)-I+A)
    cc.c_ZZ = dalek_G2::mul_by_a((T1-X1)*(T2+X2)-I+A); // dalek_param_twist_coeff_a is 1*X for us

    cc.c_XY = X1*Z2-X2*Z1+F;          // c_XY = X1*Z2-X2*Z1+F
    cc.c_XZ = (Y1-T1)*(Y2+T2)-B+I-H;  // c_XZ = (Y1-T1)*(Y2+T2)-B+I-H
    current.X = E*F;                  // X3   = E*F
    current.Y = G*H;                  // Y3   = G*H
    current.Z = F*G;                  // Z3   = F*G
    current.T = E*H;                  // T3   = E*H

#ifdef DEBUG
    current.test_invariant();
#endif
}

void mixed_addition_step_for_flipped_miller_loop(const extended_dalek_G2_projective &base,
                                                 extended_dalek_G2_projective &current,
                                                 dalek_Fq3_conic_coefficients &cc)
{
    const dalek_Fq3 &X1 = current.X, &Y1 = current.Y, &Z1 = current.Z, &T1 = current.T;
    const dalek_Fq3 &X2 = base.X, &Y2 =  base.Y, &T2 = base.T;

    const dalek_Fq3 A = X1*X2;               // A    = X1*X2
    const dalek_Fq3 B = Y1*Y2;               // B    = Y1*Y2
    const dalek_Fq3 C = Z1*T2;               // C    = Z1*T2
    const dalek_Fq3 E = T1+C;                // E    = T1+C
    const dalek_Fq3 F = (X1-Y1)*(X2+Y2)+B-A; // F    = (X1-Y1)*(X2+Y2)+B-A
    // G = B + twisted_dalek_a * A
    const dalek_Fq3 G = B + dalek_G2::mul_by_a(A); // dalek_param_twist_coeff_a is 1*X for us
    const dalek_Fq3 H = T1-C;                // H    = T1-C
    const dalek_Fq3 I = T1*T2;               // I    = T1*T2

    // c_ZZ = delta_3* ((T1-X1)*(T2+X2)-I+A)
    cc.c_ZZ = dalek_G2::mul_by_a((T1-X1)*(T2+X2)-I+A); // dalek_param_twist_coeff_a is 1*X for us

    cc.c_XY = X1-X2*Z1+F;             // c_XY = X1*Z2-X2*Z1+F
    cc.c_XZ = (Y1-T1)*(Y2+T2)-B+I-H;  // c_XZ = (Y1-T1)*(Y2+T2)-B+I-H
    current.X = E*F;                  // X3   = E*F
    current.Y = G*H;                  // Y3   = G*H
    current.Z = F*G;                  // Z3   = F*G
    current.T = E*H;                  // T3   = E*H

#ifdef DEBUG
    current.test_invariant();
#endif
}

dalek_ate_G1_precomp dalek_ate_precompute_G1(const dalek_G1& P)
{
    enter_block("Call to dalek_ate_precompute_G1");
    dalek_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();
    dalek_ate_G1_precomp result;
    result.P_XY = Pcopy.X*Pcopy.Y;
    result.P_XZ = Pcopy.X; // P.X * P.Z but P.Z = 1
    result.P_ZZplusYZ = (dalek_Fq::one() + Pcopy.Y); // (P.Z + P.Y) * P.Z but P.Z = 1
    leave_block("Call to dalek_ate_precompute_G1");
    return result;
}

dalek_ate_G2_precomp dalek_ate_precompute_G2(const dalek_G2& Q)
{
    enter_block("Call to dalek_ate_precompute_G2");
    const bigint<dalek_Fr::num_limbs> &loop_count = dalek_ate_loop_count;
    dalek_ate_G2_precomp result;

    dalek_G2 Qcopy(Q);
    Qcopy.to_affine_coordinates();

    extended_dalek_G2_projective Q_ext;
    Q_ext.X = Qcopy.X;
    Q_ext.Y = Qcopy.Y;
    Q_ext.Z = Qcopy.Z;
    Q_ext.T = Qcopy.X*Qcopy.Y;

    extended_dalek_G2_projective R = Q_ext;

    bool found_one = false;
    for (long i = loop_count.max_bits()-1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        dalek_Fq3_conic_coefficients cc;
        doubling_step_for_flipped_miller_loop(R, cc);
        result.push_back(cc);
        if (bit)
        {
            mixed_addition_step_for_flipped_miller_loop(Q_ext, R, cc);
            result.push_back(cc);
        }
    }

    leave_block("Call to dalek_ate_precompute_G2");
    return result;
}

dalek_Fq6 dalek_ate_miller_loop(const dalek_ate_G1_precomp &prec_P,
                                    const dalek_ate_G2_precomp &prec_Q)
{
    enter_block("Call to dalek_ate_miller_loop");
    const bigint<dalek_Fr::num_limbs> &loop_count = dalek_ate_loop_count;

    dalek_Fq6 f = dalek_Fq6::one();

    bool found_one = false;
    size_t idx = 0;
    for (long i = loop_count.max_bits()-1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           dalek_param_p (skipping leading zeros) in MSB to LSB
           order */
        dalek_Fq3_conic_coefficients cc = prec_Q[idx++];

        dalek_Fq6 g_RR_at_P = dalek_Fq6(prec_P.P_XY * cc.c_XY + prec_P.P_XZ * cc.c_XZ,
                                            prec_P.P_ZZplusYZ * cc.c_ZZ);
        f = f.squared() * g_RR_at_P;
        if (bit)
        {
            cc = prec_Q[idx++];
            dalek_Fq6 g_RQ_at_P = dalek_Fq6(prec_P.P_ZZplusYZ * cc.c_ZZ,
                                                prec_P.P_XY * cc.c_XY + prec_P.P_XZ * cc.c_XZ);
            f = f * g_RQ_at_P;
        }
    }
    leave_block("Call to dalek_ate_miller_loop");

    return f;
}

dalek_Fq6 dalek_ate_double_miller_loop(const dalek_ate_G1_precomp &prec_P1,
                                           const dalek_ate_G2_precomp &prec_Q1,
                                           const dalek_ate_G1_precomp &prec_P2,
                                           const dalek_ate_G2_precomp &prec_Q2)
{
    enter_block("Call to dalek_ate_double_miller_loop");
    const bigint<dalek_Fr::num_limbs> &loop_count = dalek_ate_loop_count;

    dalek_Fq6 f = dalek_Fq6::one();

    bool found_one = false;
    size_t idx = 0;
    for (long i = loop_count.max_bits()-1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           dalek_param_p (skipping leading zeros) in MSB to LSB
           order */
        dalek_Fq3_conic_coefficients cc1 = prec_Q1[idx];
        dalek_Fq3_conic_coefficients cc2 = prec_Q2[idx];
        ++idx;

        dalek_Fq6 g_RR_at_P1 = dalek_Fq6(prec_P1.P_XY * cc1.c_XY + prec_P1.P_XZ * cc1.c_XZ,
                                             prec_P1.P_ZZplusYZ * cc1.c_ZZ);

        dalek_Fq6 g_RR_at_P2 = dalek_Fq6(prec_P2.P_XY * cc2.c_XY + prec_P2.P_XZ * cc2.c_XZ,
                                             prec_P2.P_ZZplusYZ * cc2.c_ZZ);
        f = f.squared() * g_RR_at_P1 * g_RR_at_P2;

        if (bit)
        {
            cc1 = prec_Q1[idx];
            cc2 = prec_Q2[idx];
            ++idx;
            dalek_Fq6 g_RQ_at_P1 = dalek_Fq6(prec_P1.P_ZZplusYZ * cc1.c_ZZ,
                                                 prec_P1.P_XY * cc1.c_XY + prec_P1.P_XZ * cc1.c_XZ);
            dalek_Fq6 g_RQ_at_P2 = dalek_Fq6(prec_P2.P_ZZplusYZ * cc2.c_ZZ,
                                                 prec_P2.P_XY * cc2.c_XY + prec_P2.P_XZ * cc2.c_XZ);
            f = f * g_RQ_at_P1 * g_RQ_at_P2;
        }
    }
    leave_block("Call to dalek_ate_double_miller_loop");

    return f;
}

dalek_Fq6 dalek_ate_pairing(const dalek_G1& P, const dalek_G2 &Q)
{
    enter_block("Call to dalek_ate_pairing");
    dalek_ate_G1_precomp prec_P = dalek_ate_precompute_G1(P);
    dalek_ate_G2_precomp prec_Q = dalek_ate_precompute_G2(Q);
    dalek_Fq6 result = dalek_ate_miller_loop(prec_P, prec_Q);
    leave_block("Call to dalek_ate_pairing");
    return result;
}

dalek_GT dalek_ate_reduced_pairing(const dalek_G1 &P, const dalek_G2 &Q)
{
    enter_block("Call to dalek_ate_reduced_pairing");
    const dalek_Fq6 f = dalek_ate_pairing(P, Q);
    const dalek_GT result = dalek_final_exponentiation(f);
    leave_block("Call to dalek_ate_reduced_pairing");
    return result;
}

dalek_G1_precomp dalek_precompute_G1(const dalek_G1& P)
{
    return dalek_ate_precompute_G1(P);
}

dalek_G2_precomp dalek_precompute_G2(const dalek_G2& Q)
{
    return dalek_ate_precompute_G2(Q);
}

dalek_Fq6 dalek_miller_loop(const dalek_G1_precomp &prec_P,
                                const dalek_G2_precomp &prec_Q)
{
    return dalek_ate_miller_loop(prec_P, prec_Q);
}

dalek_Fq6 dalek_double_miller_loop(const dalek_G1_precomp &prec_P1,
                                       const dalek_G2_precomp &prec_Q1,
                                       const dalek_G1_precomp &prec_P2,
                                       const dalek_G2_precomp &prec_Q2)
{
    return dalek_ate_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

dalek_Fq6 dalek_pairing(const dalek_G1& P,
                            const dalek_G2 &Q)
{
    return dalek_ate_pairing(P, Q);
}

dalek_GT dalek_reduced_pairing(const dalek_G1 &P,
                                   const dalek_G2 &Q)
{
    return dalek_ate_reduced_pairing(P, Q);
}
} // libff
