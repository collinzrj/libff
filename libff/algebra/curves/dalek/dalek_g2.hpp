/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef DALEK_G2_HPP_
#define DALEK_G2_HPP_
#include <iostream>
#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/dalek/dalek_init.hpp>

namespace libff {

class dalek_G2;
std::ostream& operator<<(std::ostream &, const dalek_G2&);
std::istream& operator>>(std::istream &, dalek_G2&);

class dalek_G2 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;

    static dalek_G2 G2_zero;
    static dalek_G2 G2_one;

    dalek_Fq3 X, Y, Z;
    dalek_G2();
private:
    dalek_G2(const dalek_Fq3& X, const dalek_Fq3& Y, const dalek_Fq3& Z) : X(X), Y(Y), Z(Z) {};
public:
    static dalek_Fq3 mul_by_a(const dalek_Fq3 &elt);
    static dalek_Fq3 mul_by_d(const dalek_Fq3 &elt);
    typedef dalek_Fq base_field;
    typedef dalek_Fq3 twist_field;
    typedef dalek_Fr scalar_field;

    // using inverted coordinates
    dalek_G2(const dalek_Fq3& X, const dalek_Fq3& Y) : X(Y), Y(X), Z(X*Y) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const dalek_G2 &other) const;
    bool operator!=(const dalek_G2 &other) const;

    dalek_G2 operator+(const dalek_G2 &other) const;
    dalek_G2 operator-() const;
    dalek_G2 operator-(const dalek_G2 &other) const;

    dalek_G2 add(const dalek_G2 &other) const;
    dalek_G2 mixed_add(const dalek_G2 &other) const;
    dalek_G2 dbl() const;
    dalek_G2 mul_by_q() const;

    bool is_well_formed() const;

    static dalek_G2 zero();
    static dalek_G2 one();
    static dalek_G2 random_element();

    static size_t size_in_bits() { return twist_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const dalek_G2 &g);
    friend std::istream& operator>>(std::istream &in, dalek_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<dalek_G2> &vec);
};

template<mp_size_t m>
dalek_G2 operator*(const bigint<m> &lhs, const dalek_G2 &rhs)
{
    return scalar_mul<dalek_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
dalek_G2 operator*(const Fp_model<m, modulus_p> &lhs, const dalek_G2 &rhs)
{
   return scalar_mul<dalek_G2, m>(rhs, lhs.as_bigint());
}

} // libff
#endif // DALEK_G2_HPP_
