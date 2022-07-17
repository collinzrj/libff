/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef DALEK_G1_HPP_
#define DALEK_G1_HPP_
#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/dalek/dalek_init.hpp>

namespace libff {

class dalek_G1;
std::ostream& operator<<(std::ostream &, const dalek_G1&);
std::istream& operator>>(std::istream &, dalek_G1&);

class dalek_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static dalek_G1 G1_zero;
    static dalek_G1 G1_one;

    dalek_Fq X, Y, Z;
    dalek_G1();
private:
    dalek_G1(const dalek_Fq& X, const dalek_Fq& Y, const dalek_Fq& Z) : X(X), Y(Y), Z(Z) {};

public:
    typedef dalek_Fq base_field;
    typedef dalek_Fr scalar_field;
    // using inverted coordinates
    dalek_G1(const dalek_Fq& X, const dalek_Fq& Y) : X(Y), Y(X), Z(X*Y) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const dalek_G1 &other) const;
    bool operator!=(const dalek_G1 &other) const;

    dalek_G1 operator+(const dalek_G1 &other) const;
    dalek_G1 operator-() const;
    dalek_G1 operator-(const dalek_G1 &other) const;

    dalek_G1 add(const dalek_G1 &other) const;
    dalek_G1 mixed_add(const dalek_G1 &other) const;
    dalek_G1 dbl() const;

    bool is_well_formed() const;

    static dalek_G1 zero();
    static dalek_G1 one();
    static dalek_G1 random_element();

    static size_t size_in_bits() { return dalek_Fq::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const dalek_G1 &g);
    friend std::istream& operator>>(std::istream &in, dalek_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<dalek_G1> &vec);
};

template<mp_size_t m>
dalek_G1 operator*(const bigint<m> &lhs, const dalek_G1 &rhs)
{
    return scalar_mul<dalek_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
dalek_G1 operator*(const Fp_model<m,modulus_p> &lhs, const dalek_G1 &rhs)
{
    return scalar_mul<dalek_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<dalek_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<dalek_G1> &v);

} // libff
#endif // DALEK_G1_HPP_
