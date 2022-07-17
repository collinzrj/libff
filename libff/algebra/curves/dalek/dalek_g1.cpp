/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/dalek/dalek_g1.hpp>

namespace libff {

#ifdef PROFILE_OP_COUNTS
long long dalek_G1::add_cnt = 0;
long long dalek_G1::dbl_cnt = 0;
#endif

std::vector<size_t> dalek_G1::wnaf_window_table;
std::vector<size_t> dalek_G1::fixed_base_exp_window_table;
dalek_G1 dalek_G1::G1_zero;
dalek_G1 dalek_G1::G1_one;

dalek_G1::dalek_G1()
{
    this->X = G1_zero.X;
    this->Y = G1_zero.Y;
    this->Z = G1_zero.Z;
}

void dalek_G1::print() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        dalek_G1 copy(*this);
        copy.to_affine_coordinates();
        gmp_printf("(%Nd , %Nd)\n",
                   copy.X.as_bigint().data, dalek_Fq::num_limbs,
                   copy.Y.as_bigint().data, dalek_Fq::num_limbs);
    }
}

void dalek_G1::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        gmp_printf("(%Nd : %Nd : %Nd)\n",
                   this->X.as_bigint().data, dalek_Fq::num_limbs,
                   this->Y.as_bigint().data, dalek_Fq::num_limbs,
                   this->Z.as_bigint().data, dalek_Fq::num_limbs);
    }
}

void dalek_G1::to_affine_coordinates()
{
    if (this->is_zero())
    {
        this->X = dalek_Fq::zero();
        this->Y = dalek_Fq::one();
        this->Z = dalek_Fq::one();
    }
    else
    {
        // go from inverted coordinates to projective coordinates
        dalek_Fq tX = this->Y * this->Z;
        dalek_Fq tY = this->X * this->Z;
        dalek_Fq tZ = this->X * this->Y;
        // go from projective coordinates to affine coordinates
        dalek_Fq tZ_inv = tZ.inverse();
        this->X = tX * tZ_inv;
        this->Y = tY * tZ_inv;
        this->Z = dalek_Fq::one();
    }
}

void dalek_G1::to_special()
{
    if (this->Z.is_zero())
    {
        return;
    }

#ifdef DEBUG
    const dalek_G1 copy(*this);
#endif

    dalek_Fq Z_inv = this->Z.inverse();
    this->X = this->X * Z_inv;
    this->Y = this->Y * Z_inv;
    this->Z = dalek_Fq::one();

#ifdef DEBUG
    assert((*this) == copy);
#endif
}

bool dalek_G1::is_special() const
{
    return (this->is_zero() || this->Z == dalek_Fq::one());
}

bool dalek_G1::is_zero() const
{
    return (this->Y.is_zero() && this->Z.is_zero());
}

bool dalek_G1::operator==(const dalek_G1 &other) const
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

bool dalek_G1::operator!=(const dalek_G1& other) const
{
    return !(operator==(other));
}

dalek_G1 dalek_G1::operator+(const dalek_G1 &other) const
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

dalek_G1 dalek_G1::operator-() const
{
    return dalek_G1(-(this->X), this->Y, this->Z);
}


dalek_G1 dalek_G1::operator-(const dalek_G1 &other) const
{
    return (*this) + (-other);
}

dalek_G1 dalek_G1::add(const dalek_G1 &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-dalek-inverted.html#addition-add-2007-bl

    dalek_Fq A = (this->Z) * (other.Z);                   // A = Z1*Z2
    dalek_Fq B = dalek_coeff_d * A.squared();           // B = d*A^2
    dalek_Fq C = (this->X) * (other.X);                   // C = X1*X2
    dalek_Fq D = (this->Y) * (other.Y);                   // D = Y1*Y2
    dalek_Fq E = C * D;                                   // E = C*D
    dalek_Fq H = C - D;                                   // H = C-D
    dalek_Fq I = (this->X+this->Y)*(other.X+other.Y)-C-D; // I = (X1+Y1)*(X2+Y2)-C-D
    dalek_Fq X3 = (E+B)*H;                                // X3 = c*(E+B)*H
    dalek_Fq Y3 = (E-B)*I;                                // Y3 = c*(E-B)*I
    dalek_Fq Z3 = A*H*I;                                  // Z3 = A*H*I

    return dalek_G1(X3, Y3, Z3);
}

dalek_G1 dalek_G1::mixed_add(const dalek_G1 &other) const
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

    dalek_Fq A = this->Z;                                 // A = Z1
    dalek_Fq B = dalek_coeff_d * A.squared();           // B = d*A^2
    dalek_Fq C = (this->X) * (other.X);                   // C = X1*X2
    dalek_Fq D = (this->Y) * (other.Y);                   // D = Y1*Y2
    dalek_Fq E = C * D;                                   // E = C*D
    dalek_Fq H = C - D;                                   // H = C-D
    dalek_Fq I = (this->X+this->Y)*(other.X+other.Y)-C-D; // I = (X1+Y1)*(X2+Y2)-C-D
    dalek_Fq X3 = (E+B)*H;                                // X3 = c*(E+B)*H
    dalek_Fq Y3 = (E-B)*I;                                // Y3 = c*(E-B)*I
    dalek_Fq Z3 = A*H*I;                                  // Z3 = A*H*I

    return dalek_G1(X3, Y3, Z3);
}

dalek_G1 dalek_G1::dbl() const
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
        // http://www.hyperelliptic.org/EFD/g1p/auto-dalek-inverted.html#doubling-dbl-2007-bl

        dalek_Fq A = (this->X).squared();                      // A = X1^2
        dalek_Fq B = (this->Y).squared();                      // B = Y1^2
        dalek_Fq C = A+B;                                      // C = A+B
        dalek_Fq D = A-B;                                      // D = A-B
        dalek_Fq E = (this->X+this->Y).squared()-C;            // E = (X1+Y1)^2-C
        dalek_Fq X3 = C*D;                                     // X3 = C*D
        dalek_Fq dZZ = dalek_coeff_d * this->Z.squared();
        dalek_Fq Y3 = E*(C-dZZ-dZZ);                           // Y3 = E*(C-2*d*Z1^2)
        dalek_Fq Z3 = D*E;                                     // Z3 = D*E

        return dalek_G1(X3, Y3, Z3);
    }
}

bool dalek_G1::is_well_formed() const
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
        dalek_Fq X2 = this->X.squared();
        dalek_Fq Y2 = this->Y.squared();
        dalek_Fq Z2 = this->Z.squared();

        // for G1 a = 1
        return (Z2 * (Y2 + X2 - dalek_coeff_d * Z2) == X2 * Y2);
    }
}

dalek_G1 dalek_G1::zero()
{
    return G1_zero;
}

dalek_G1 dalek_G1::one()
{
    return G1_one;
}

dalek_G1 dalek_G1::random_element()
{
    return dalek_Fr::random_element().as_bigint() * G1_one;
}

std::ostream& operator<<(std::ostream &out, const dalek_G1 &g)
{
    dalek_G1 copy(g);
    copy.to_affine_coordinates();
#ifdef NO_PT_COMPRESSION
    out << copy.X << OUTPUT_SEPARATOR << copy.Y;
#else
    /* storing LSB of Y */
    out << copy.X << OUTPUT_SEPARATOR << (copy.Y.as_bigint().data[0] & 1);
#endif

    return out;
}

std::istream& operator>>(std::istream &in, dalek_G1 &g)
{
    dalek_Fq tX, tY;

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

    dalek_Fq tX2 = tX.squared();
    dalek_Fq tY2 = (dalek_Fq::one() - tX2) * // a = 1 for G1 (not a twist)
        (dalek_Fq::one() - dalek_coeff_d * tX2).inverse();
    tY = tY2.sqrt();

    if ((tY.as_bigint().data[0] & 1) != Y_lsb)
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

std::ostream& operator<<(std::ostream& out, const std::vector<dalek_G1> &v)
{
    out << v.size() << "\n";
    for (const dalek_G1& t : v)
    {
        out << t << OUTPUT_NEWLINE;
    }

    return out;
}

std::istream& operator>>(std::istream& in, std::vector<dalek_G1> &v)
{
    v.clear();

    size_t s;
    in >> s;
    v.reserve(s);
    consume_newline(in);

    for (size_t i = 0; i < s; ++i)
    {
        dalek_G1 g;
        in >> g;
        v.emplace_back(g);
        consume_OUTPUT_NEWLINE(in);
    }

    return in;
}

void dalek_G1::batch_to_special_all_non_zeros(std::vector<dalek_G1> &vec)
{
    std::vector<dalek_Fq> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z);
    }
    batch_invert<dalek_Fq>(Z_vec);

    const dalek_Fq one = dalek_Fq::one();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        vec[i].X = vec[i].X * Z_vec[i];
        vec[i].Y = vec[i].Y * Z_vec[i];
        vec[i].Z = one;
    }
}

} // libff
