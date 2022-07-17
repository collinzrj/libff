/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/dalek/dalek_g2.hpp>

namespace libff {

#ifdef PROFILE_OP_COUNTS
long long dalek_G2::add_cnt = 0;
long long dalek_G2::dbl_cnt = 0;
#endif

std::vector<size_t> dalek_G2::wnaf_window_table;
std::vector<size_t> dalek_G2::fixed_base_exp_window_table;

dalek_G2 dalek_G2::G2_zero;
dalek_G2 dalek_G2::G2_one;

dalek_G2::dalek_G2()
{
    this->X = G2_zero.X;
    this->Y = G2_zero.Y;
    this->Z = G2_zero.Z;
}

dalek_Fq3 dalek_G2::mul_by_a(const dalek_Fq3 &elt)
{
	// should be
	//  dalek_Fq3(dalek_twist_mul_by_a_c0 * elt.c2, dalek_twist_mul_by_a_c1 * elt.c0, dalek_twist_mul_by_a_c2 * elt.c1)
	// but optimizing the fact that dalek_twist_mul_by_a_c1 = dalek_twist_mul_by_a_c2 = 1
    return dalek_Fq3(dalek_twist_mul_by_a_c0 * elt.c2, elt.c0, elt.c1);
}

dalek_Fq3 dalek_G2::mul_by_d(const dalek_Fq3 &elt)
{
	return dalek_Fq3(dalek_twist_mul_by_d_c0 * elt.c2, dalek_twist_mul_by_d_c1 * elt.c0, dalek_twist_mul_by_d_c2 * elt.c1);
}

void dalek_G2::print() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        dalek_G2 copy(*this);
        copy.to_affine_coordinates();
        gmp_printf("(%Nd*z^2 + %Nd*z + %Nd , %Nd*z^2 + %Nd*z + %Nd)\n",
                   copy.X.c2.as_bigint().data, dalek_Fq::num_limbs,
                   copy.X.c1.as_bigint().data, dalek_Fq::num_limbs,
                   copy.X.c0.as_bigint().data, dalek_Fq::num_limbs,
                   copy.Y.c2.as_bigint().data, dalek_Fq::num_limbs,
                   copy.Y.c1.as_bigint().data, dalek_Fq::num_limbs,
                   copy.Y.c0.as_bigint().data, dalek_Fq::num_limbs);
    }
}

void dalek_G2::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        gmp_printf("(%Nd*z^2 + %Nd*z + %Nd : %Nd*z^2 + %Nd*z + %Nd : %Nd*z^2 + %Nd*z + %Nd)\n",
                   this->X.c2.as_bigint().data, dalek_Fq::num_limbs,
                   this->X.c1.as_bigint().data, dalek_Fq::num_limbs,
                   this->X.c0.as_bigint().data, dalek_Fq::num_limbs,
                   this->Y.c2.as_bigint().data, dalek_Fq::num_limbs,
                   this->Y.c1.as_bigint().data, dalek_Fq::num_limbs,
                   this->Y.c0.as_bigint().data, dalek_Fq::num_limbs,
                   this->Z.c2.as_bigint().data, dalek_Fq::num_limbs,
                   this->Z.c1.as_bigint().data, dalek_Fq::num_limbs,
                   this->Z.c0.as_bigint().data, dalek_Fq::num_limbs);
    }
}

void dalek_G2::to_affine_coordinates()
{
    if (this->is_zero())
    {
        this->X = dalek_Fq3::zero();
        this->Y = dalek_Fq3::one();
        this->Z = dalek_Fq3::one();
    }
    else
    {
        // go from inverted coordinates to projective coordinates
        dalek_Fq3 tX = this->Y * this->Z;
        dalek_Fq3 tY = this->X * this->Z;
        dalek_Fq3 tZ = this->X * this->Y;
        // go from projective coordinates to affine coordinates
        dalek_Fq3 tZ_inv = tZ.inverse();
        this->X = tX * tZ_inv;
        this->Y = tY * tZ_inv;
        this->Z = dalek_Fq3::one();
    }
}

void dalek_G2::to_special()
{
    if (this->Z.is_zero())
    {
        return;
    }

#ifdef DEBUG
    const dalek_G2 copy(*this);
#endif

    dalek_Fq3 Z_inv = this->Z.inverse();
    this->X = this->X * Z_inv;
    this->Y = this->Y * Z_inv;
    this->Z = dalek_Fq3::one();

#ifdef DEBUG
    assert((*this) == copy);
#endif
}

bool dalek_G2::is_special() const
{
    return (this->is_zero() || this->Z == dalek_Fq3::one());
}

bool dalek_G2::is_zero() const
{
    return (this->Y.is_zero() && this->Z.is_zero());
}

bool dalek_G2::operator==(const dalek_G2 &other) const
{
    if (this->is_zero())
    {
        return other.is_zero();
    }

    if (other.is_zero())
    {
        return false;
    }

    /* now neither is O */

    // X1/Z1 = X2/Z2 <=> X1*Z2 = X2*Z1
    if ((this->X * other.Z) != (other.X * this->Z))
    {
        return false;
    }

    // Y1/Z1 = Y2/Z2 <=> Y1*Z2 = Y2*Z1
    if ((this->Y * other.Z) != (other.Y * this->Z))
    {
        return false;
    }

    return true;
}

bool dalek_G2::operator!=(const dalek_G2& other) const
{
    return !(operator==(other));
}

dalek_G2 dalek_G2::operator+(const dalek_G2 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return (*this);
    }

    return this->add(other);
}

dalek_G2 dalek_G2::operator-() const
{
    return dalek_G2(-(this->X), this->Y, this->Z);
}


dalek_G2 dalek_G2::operator-(const dalek_G2 &other) const
{
    return (*this) + (-other);
}

dalek_G2 dalek_G2::add(const dalek_G2 &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-twisted-inverted.html#addition-add-2008-bbjlp

    const dalek_Fq3 A = (this->Z) * (other.Z);                       // A = Z1*Z2
    const dalek_Fq3 B = dalek_G2::mul_by_d(A.squared());           // B = d*A^2
    const dalek_Fq3 C = (this->X) * (other.X);                       // C = X1*X2
    const dalek_Fq3 D = (this->Y) * (other.Y);                       // D = Y1*Y2
    const dalek_Fq3 E = C*D;                                         // E = C*D
    const dalek_Fq3 H = C - dalek_G2::mul_by_a(D);                 // H = C-a*D
    const dalek_Fq3 I = (this->X+this->Y)*(other.X+other.Y)-C-D;     // I = (X1+Y1)*(X2+Y2)-C-D
    const dalek_Fq3 X3 = (E+B)*H;                                    // X3 = (E+B)*H
    const dalek_Fq3 Y3 = (E-B)*I;                                    // Y3 = (E-B)*I
    const dalek_Fq3 Z3 = A*H*I;                                      // Z3 = A*H*I

    return dalek_G2(X3, Y3, Z3);
}

dalek_G2 dalek_G2::mixed_add(const dalek_G2 &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

#ifdef DEBUG
    assert(other.is_special());
#endif

    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-dalek-inverted.html#addition-madd-2007-lb

    const dalek_Fq3 A = this->Z;                                     // A = Z1*Z2
    const dalek_Fq3 B = dalek_G2::mul_by_d(A.squared());           // B = d*A^2
    const dalek_Fq3 C = (this->X) * (other.X);                       // C = X1*X2
    const dalek_Fq3 D = (this->Y) * (other.Y);                       // D = Y1*Y2
    const dalek_Fq3 E = C*D;                                         // E = C*D
    const dalek_Fq3 H = C - dalek_G2::mul_by_a(D);                 // H = C-a*D
    const dalek_Fq3 I = (this->X+this->Y)*(other.X+other.Y)-C-D;     // I = (X1+Y1)*(X2+Y2)-C-D
    const dalek_Fq3 X3 = (E+B)*H;                                    // X3 = (E+B)*H
    const dalek_Fq3 Y3 = (E-B)*I;                                    // Y3 = (E-B)*I
    const dalek_Fq3 Z3 = A*H*I;                                      // Z3 = A*H*I

    return dalek_G2(X3, Y3, Z3);
}

dalek_G2 dalek_G2::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif
    if (this->is_zero())
    {
        return (*this);
    }
    else
    {
        // NOTE: does not handle O and pts of order 2,4
        // http://www.hyperelliptic.org/EFD/g1p/auto-twisted-inverted.html#doubling-dbl-2008-bbjlp

        const dalek_Fq3 A = (this->X).squared();                      // A = X1^2
        const dalek_Fq3 B = (this->Y).squared();                      // B = Y1^2
        const dalek_Fq3 U = dalek_G2::mul_by_a(B);                  // U = a*B
        const dalek_Fq3 C = A+U;                                      // C = A+U
        const dalek_Fq3 D = A-U;                                      // D = A-U
        const dalek_Fq3 E = (this->X+this->Y).squared()-A-B;          // E = (X1+Y1)^2-A-B
        const dalek_Fq3 X3 = C*D;                                     // X3 = C*D
        const dalek_Fq3 dZZ = dalek_G2::mul_by_d(this->Z.squared());
        const dalek_Fq3 Y3 = E*(C-dZZ-dZZ);                           // Y3 = E*(C-2*d*Z1^2)
        const dalek_Fq3 Z3 = D*E;                                     // Z3 = D*E

        return dalek_G2(X3, Y3, Z3);
    }
}

dalek_G2 dalek_G2::mul_by_q() const
{
    return dalek_G2((this->X).Frobenius_map(1),
                      dalek_twist_mul_by_q_Y * (this->Y).Frobenius_map(1),
                      dalek_twist_mul_by_q_Z * (this->Z).Frobenius_map(1));
}

bool dalek_G2::is_well_formed() const
{
    /* Note that point at infinity is the only special case we must check as
       inverted representation does no cover points (0, +-c) and (+-c, 0). */
    if (this->is_zero())
    {
        return true;
    }
    else
    {
        /*
          a x^2 + y^2 = 1 + d x^2 y^2

          We are using inverted, so equation we need to check is actually

          a (z/x)^2 + (z/y)^2 = 1 + d z^4 / (x^2 * y^2)
          z^2 (a y^2 + x^2 - dz^2) = x^2 y^2
        */
        dalek_Fq3 X2 = this->X.squared();
        dalek_Fq3 Y2 = this->Y.squared();
        dalek_Fq3 Z2 = this->Z.squared();
        dalek_Fq3 aY2 = dalek_G2::mul_by_a(Y2);
        dalek_Fq3 dZ2 = dalek_G2::mul_by_d(Z2);
        return (Z2 * (aY2 + X2 - dZ2) == X2 * Y2);
    }
}

dalek_G2 dalek_G2::zero()
{
    return G2_zero;
}

dalek_G2 dalek_G2::one()
{
    return G2_one;
}

dalek_G2 dalek_G2::random_element()
{
    return dalek_Fr::random_element().as_bigint() * G2_one;
}

std::ostream& operator<<(std::ostream &out, const dalek_G2 &g)
{
    dalek_G2 copy(g);
    copy.to_affine_coordinates();
#ifdef NO_PT_COMPRESSION
    out << copy.X << OUTPUT_SEPARATOR << copy.Y;
#else
    /* storing LSB of Y */
    out << copy.X << OUTPUT_SEPARATOR << (copy.Y.c0.as_bigint().data[0] & 1);
#endif
    return out;
}

std::istream& operator>>(std::istream &in, dalek_G2 &g)
{
    dalek_Fq3 tX, tY;

#ifdef NO_PT_COMPRESSION
    in >> tX;
    consume_OUTPUT_SEPARATOR(in);
    in >> tY;
#else
    /*
      a x^2 + y^2 = 1 + d x^2 y^2
      y = sqrt((1-ax^2)/(1-dx^2))
    */
    unsigned char Y_lsb;
    in >> tX;
    consume_OUTPUT_SEPARATOR(in);

    in.read((char*)&Y_lsb, 1);
    Y_lsb -= '0';

    dalek_Fq3 tX2 = tX.squared();
    const dalek_Fq3 tY2 =
        (dalek_Fq3::one() - dalek_G2::mul_by_a(tX2)) *
        (dalek_Fq3::one() - dalek_G2::mul_by_d(tX2)).inverse();
    tY = tY2.sqrt();

    if ((tY.c0.as_bigint().data[0] & 1) != Y_lsb)
    {
        tY = -tY;
    }
#endif

    // using inverted coordinates
    g.X = tY;
    g.Y = tX;
    g.Z = tX * tY;

#ifdef USE_MIXED_ADDITION
    g.to_special();
#endif

    return in;
}

void dalek_G2::batch_to_special_all_non_zeros(std::vector<dalek_G2> &vec)
{
    std::vector<dalek_Fq3> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z);
    }
    batch_invert<dalek_Fq3>(Z_vec);

    const dalek_Fq3 one = dalek_Fq3::one();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        vec[i].X = vec[i].X * Z_vec[i];
        vec[i].Y = vec[i].Y * Z_vec[i];
        vec[i].Z = one;
    }
}

} // libff
