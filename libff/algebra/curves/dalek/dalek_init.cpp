/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/dalek/dalek_g1.hpp>
#include <libff/algebra/curves/dalek/dalek_g2.hpp>
#include <libff/algebra/curves/dalek/dalek_init.hpp>

namespace libff {

bigint<dalek_r_limbs> dalek_modulus_r;
bigint<dalek_q_limbs> dalek_modulus_q;

dalek_Fq dalek_coeff_a;
dalek_Fq dalek_coeff_d;
dalek_Fq3 dalek_twist;
dalek_Fq3 dalek_twist_coeff_a;
dalek_Fq3 dalek_twist_coeff_d;
dalek_Fq dalek_twist_mul_by_a_c0;
dalek_Fq dalek_twist_mul_by_a_c1;
dalek_Fq dalek_twist_mul_by_a_c2;
dalek_Fq dalek_twist_mul_by_d_c0;
dalek_Fq dalek_twist_mul_by_d_c1;
dalek_Fq dalek_twist_mul_by_d_c2;
dalek_Fq dalek_twist_mul_by_q_Y;
dalek_Fq dalek_twist_mul_by_q_Z;

bigint<dalek_q_limbs> dalek_ate_loop_count;
bigint<6*dalek_q_limbs> dalek_final_exponent;
bigint<dalek_q_limbs> dalek_final_exponent_last_chunk_abs_of_w0;
bool dalek_final_exponent_last_chunk_is_w0_neg;
bigint<dalek_q_limbs> dalek_final_exponent_last_chunk_w1;

void init_dalek_params()
{
    typedef bigint<dalek_r_limbs> bigint_r;
    typedef bigint<dalek_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    //dalek_modulus_r = bigint_r("21888242871839275222246405745257275088548364400416034343698204186575808495617");
    dalek_modulus_r = bigint_r("7237005577332262213973186563042994240857116359379907606001950938285454250989");
    assert(dalek_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
	    /* bn128
        dalek_Fr::Rsquared = bigint_r("944936681149208446651664254269745548490766851729442924617792859073125903783");
        dalek_Fr::Rcubed = bigint_r("5866548545943845227489894872040244720403868105578784105281690076696998248512");
        dalek_Fr::inv = 0xc2e1f593efffffff;
	*/
        dalek_Fr::Rsquared = bigint_r("1627715501170711445284395025044413883736156588369414752970002579683115011841");
        dalek_Fr::Rcubed = bigint_r("6479107319630627712989486218399433800251238082664018703558868829414071968475");
        dalek_Fr::inv = 0xd2b51da312547e1b;
    }
    if (sizeof(mp_limb_t) == 4)
    {
	    /*
        dalek_Fr::Rsquared = bigint_r("944936681149208446651664254269745548490766851729442924617792859073125903783");
        dalek_Fr::Rcubed = bigint_r("5866548545943845227489894872040244720403868105578784105281690076696998248512");
        dalek_Fr::inv = 0xefffffff;
	*/
        dalek_Fr::Rsquared = bigint_r("1627715501170711445284395025044413883736156588369414752970002579683115011841");
        dalek_Fr::Rcubed = bigint_r("6479107319630627712989486218399433800251238082664018703558868829414071968475");
        dalek_Fr::inv = 0x12547e1b;
    }
    dalek_Fr::num_bits = 254;
    dalek_Fr::euler = bigint_r("776255515051215125618400780672310996630960448785612800");
    dalek_Fr::s = 31;
    dalek_Fr::t = bigint_r("722944284836962004768104088187507350585386575");
    dalek_Fr::t_minus_1_over_2 = bigint_r("361472142418481002384052044093753675292693287");
    dalek_Fr::multiplicative_generator = dalek_Fr("19");
    dalek_Fr::root_of_unity = dalek_Fr("695314865466598274460565335217615316274564719601897184");
    dalek_Fr::nqr = dalek_Fr("11");
    dalek_Fr::nqr_to_t = dalek_Fr("1326707053668679463752768729767248251415639579872144553");

    /* parameters for base field Fq */

    dalek_modulus_q = bigint_q("6210044120409721004947206240885978274523751269793792001");
    assert(dalek_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        dalek_Fq::Rsquared = bigint_q("5943559676554581037560514598978484097352477055348195432");
        dalek_Fq::Rcubed = bigint_q("1081560488703514202058739223469726982199727506489234349");
        dalek_Fq::inv = 0x76eb690b7fffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        dalek_Fq::Rsquared = bigint_q("5943559676554581037560514598978484097352477055348195432");
        dalek_Fq::Rcubed = bigint_q("1081560488703514202058739223469726982199727506489234349");
        dalek_Fq::inv = 0x7fffffff;
    }
    dalek_Fq::num_bits = 183;
    dalek_Fq::euler = bigint_q("3105022060204860502473603120442989137261875634896896000");
    dalek_Fq::s = 31;
    dalek_Fq::t = bigint_q("2891777139347848019072416350658041552884388375");
    dalek_Fq::t_minus_1_over_2 = bigint_q("1445888569673924009536208175329020776442194187");
    dalek_Fq::multiplicative_generator = dalek_Fq("61");
    dalek_Fq::root_of_unity = dalek_Fq("4692813029219384139894873043933463717810008194158530536");
    dalek_Fq::nqr = dalek_Fq("23");
    dalek_Fq::nqr_to_t = dalek_Fq("2626736066325740702418554487368721595489070118548299138");

    /* parameters for twist field Fq3 */

    dalek_Fq3::euler = bigint<3*dalek_q_limbs>("119744082713971502962992613191067836698205043373978948903839934564152994858051284658545502971203325031831647424413111161318314144765646525057914792711854057586688000");
    dalek_Fq3::s = 31;
    dalek_Fq3::t = bigint<3*dalek_q_limbs>("111520367408144756185815309352304634357062208814526860512643991563611659089151103662834971185031649686239331424621037357783237607000066456438894190557165125");
    dalek_Fq3::t_minus_1_over_2 = bigint<3*dalek_q_limbs>("55760183704072378092907654676152317178531104407263430256321995781805829544575551831417485592515824843119665712310518678891618803500033228219447095278582562");
    dalek_Fq3::non_residue = dalek_Fq("61");
    dalek_Fq3::nqr = dalek_Fq3(dalek_Fq("23"),dalek_Fq("0"),dalek_Fq("0"));
    dalek_Fq3::nqr_to_t = dalek_Fq3(dalek_Fq("104810943629412208121981114244673004633270996333237516"),dalek_Fq("0"),dalek_Fq("0"));
    dalek_Fq3::Frobenius_coeffs_c1[0] = dalek_Fq("1");
    dalek_Fq3::Frobenius_coeffs_c1[1] = dalek_Fq("1073752683758513276629212192812154536507607213288832061");
    dalek_Fq3::Frobenius_coeffs_c1[2] = dalek_Fq("5136291436651207728317994048073823738016144056504959939");
    dalek_Fq3::Frobenius_coeffs_c2[0] = dalek_Fq("1");
    dalek_Fq3::Frobenius_coeffs_c2[1] = dalek_Fq("5136291436651207728317994048073823738016144056504959939");
    dalek_Fq3::Frobenius_coeffs_c2[2] = dalek_Fq("1073752683758513276629212192812154536507607213288832061");

    /* parameters for Fq6 */

    dalek_Fq6::non_residue = dalek_Fq("61");
    dalek_Fq6::Frobenius_coeffs_c1[0] = dalek_Fq("1");
    dalek_Fq6::Frobenius_coeffs_c1[1] = dalek_Fq("1073752683758513276629212192812154536507607213288832062");
    dalek_Fq6::Frobenius_coeffs_c1[2] = dalek_Fq("1073752683758513276629212192812154536507607213288832061");
    dalek_Fq6::Frobenius_coeffs_c1[3] = dalek_Fq("6210044120409721004947206240885978274523751269793792000");
    dalek_Fq6::Frobenius_coeffs_c1[4] = dalek_Fq("5136291436651207728317994048073823738016144056504959939");
    dalek_Fq6::Frobenius_coeffs_c1[5] = dalek_Fq("5136291436651207728317994048073823738016144056504959940");
    dalek_Fq6::my_Fp2::non_residue = dalek_Fq3::non_residue;

    /* choice of Edwards curve and its twist */

    dalek_coeff_a = dalek_Fq::one();
    dalek_coeff_d = dalek_Fq("600581931845324488256649384912508268813600056237543024");
    dalek_twist = dalek_Fq3(dalek_Fq::zero(), dalek_Fq::one(), dalek_Fq::zero());
    dalek_twist_coeff_a = dalek_coeff_a * dalek_twist;
    dalek_twist_coeff_d = dalek_coeff_d * dalek_twist;
    dalek_twist_mul_by_a_c0 = dalek_coeff_a * dalek_Fq3::non_residue;
    dalek_twist_mul_by_a_c1 = dalek_coeff_a;
    dalek_twist_mul_by_a_c2 = dalek_coeff_a;
    dalek_twist_mul_by_d_c0 = dalek_coeff_d * dalek_Fq3::non_residue;
    dalek_twist_mul_by_d_c1 = dalek_coeff_d;
    dalek_twist_mul_by_d_c2 = dalek_coeff_d;
    dalek_twist_mul_by_q_Y = dalek_Fq("1073752683758513276629212192812154536507607213288832062");
    dalek_twist_mul_by_q_Z = dalek_Fq("1073752683758513276629212192812154536507607213288832062");

    /* choice of group G1 */

    dalek_G1::G1_zero = dalek_G1(dalek_Fq::zero(),
                                     dalek_Fq::one());
    dalek_G1::G1_one = dalek_G1(dalek_Fq("3713709671941291996998665608188072510389821008693530490"),
                                    dalek_Fq("4869953702976555123067178261685365085639705297852816679"));

    dalek_G1::wnaf_window_table.resize(0);
    dalek_G1::wnaf_window_table.push_back(9);
    dalek_G1::wnaf_window_table.push_back(14);
    dalek_G1::wnaf_window_table.push_back(24);
    dalek_G1::wnaf_window_table.push_back(117);

    dalek_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.10]
    dalek_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.10, 9.69]
    dalek_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [9.69, 25.21]
    dalek_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.21, 60.00]
    dalek_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.00, 149.33]
    dalek_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [149.33, 369.61]
    dalek_G1::fixed_base_exp_window_table.push_back(149);
    // window 7 is unbeaten in [369.61, 849.07]
    dalek_G1::fixed_base_exp_window_table.push_back(370);
    // window 8 is unbeaten in [849.07, 1764.94]
    dalek_G1::fixed_base_exp_window_table.push_back(849);
    // window 9 is unbeaten in [1764.94, 4429.59]
    dalek_G1::fixed_base_exp_window_table.push_back(1765);
    // window 10 is unbeaten in [4429.59, 13388.78]
    dalek_G1::fixed_base_exp_window_table.push_back(4430);
    // window 11 is unbeaten in [13388.78, 15368.00]
    dalek_G1::fixed_base_exp_window_table.push_back(13389);
    // window 12 is unbeaten in [15368.00, 74912.07]
    dalek_G1::fixed_base_exp_window_table.push_back(15368);
    // window 13 is unbeaten in [74912.07, 438107.20]
    dalek_G1::fixed_base_exp_window_table.push_back(74912);
    // window 14 is never the best
    dalek_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [438107.20, 1045626.18]
    dalek_G1::fixed_base_exp_window_table.push_back(438107);
    // window 16 is never the best
    dalek_G1::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [1045626.18, 1577434.48]
    dalek_G1::fixed_base_exp_window_table.push_back(1045626);
    // window 18 is unbeaten in [1577434.48, 17350594.23]
    dalek_G1::fixed_base_exp_window_table.push_back(1577434);
    // window 19 is never the best
    dalek_G1::fixed_base_exp_window_table.push_back(0);
    // window 20 is never the best
    dalek_G1::fixed_base_exp_window_table.push_back(0);
    // window 21 is unbeaten in [17350594.23, inf]
    dalek_G1::fixed_base_exp_window_table.push_back(17350594);
    // window 22 is never the best
    dalek_G1::fixed_base_exp_window_table.push_back(0);

    /* choice of group G2 */

    dalek_G2::G2_zero = dalek_G2(dalek_Fq3::zero(),
                                     dalek_Fq3::one());
    dalek_G2::G2_one = dalek_G2(dalek_Fq3(dalek_Fq("4531683359223370252210990718516622098304721701253228128"),
                                                dalek_Fq("5339624155305731263217400504407647531329993548123477368"),
                                                dalek_Fq("3964037981777308726208525982198654699800283729988686552")),
                                    dalek_Fq3(dalek_Fq("364634864866983740775341816274081071386963546650700569"),
                                                dalek_Fq("3264380230116139014996291397901297105159834497864380415"),
                                                dalek_Fq("3504781284999684163274269077749440837914479176282903747")));

    dalek_G2::wnaf_window_table.resize(0);
    dalek_G2::wnaf_window_table.push_back(6);
    dalek_G2::wnaf_window_table.push_back(12);
    dalek_G2::wnaf_window_table.push_back(42);
    dalek_G2::wnaf_window_table.push_back(97);

    dalek_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.74]
    dalek_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.74, 10.67]
    dalek_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.67, 25.53]
    dalek_G2::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [25.53, 60.67]
    dalek_G2::fixed_base_exp_window_table.push_back(26);
    // window 5 is unbeaten in [60.67, 145.77]
    dalek_G2::fixed_base_exp_window_table.push_back(61);
    // window 6 is unbeaten in [145.77, 356.76]
    dalek_G2::fixed_base_exp_window_table.push_back(146);
    // window 7 is unbeaten in [356.76, 823.08]
    dalek_G2::fixed_base_exp_window_table.push_back(357);
    // window 8 is unbeaten in [823.08, 1589.45]
    dalek_G2::fixed_base_exp_window_table.push_back(823);
    // window 9 is unbeaten in [1589.45, 4135.70]
    dalek_G2::fixed_base_exp_window_table.push_back(1589);
    // window 10 is unbeaten in [4135.70, 14297.74]
    dalek_G2::fixed_base_exp_window_table.push_back(4136);
    // window 11 is unbeaten in [14297.74, 16744.85]
    dalek_G2::fixed_base_exp_window_table.push_back(14298);
    // window 12 is unbeaten in [16744.85, 51768.98]
    dalek_G2::fixed_base_exp_window_table.push_back(16745);
    // window 13 is unbeaten in [51768.98, 99811.01]
    dalek_G2::fixed_base_exp_window_table.push_back(51769);
    // window 14 is unbeaten in [99811.01, 193306.72]
    dalek_G2::fixed_base_exp_window_table.push_back(99811);
    // window 15 is unbeaten in [193306.72, 907184.68]
    dalek_G2::fixed_base_exp_window_table.push_back(193307);
    // window 16 is never the best
    dalek_G2::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [907184.68, 1389682.59]
    dalek_G2::fixed_base_exp_window_table.push_back(907185);
    // window 18 is unbeaten in [1389682.59, 6752695.74]
    dalek_G2::fixed_base_exp_window_table.push_back(1389683);
    // window 19 is never the best
    dalek_G2::fixed_base_exp_window_table.push_back(0);
    // window 20 is unbeaten in [6752695.74, 193642894.51]
    dalek_G2::fixed_base_exp_window_table.push_back(6752696);
    // window 21 is unbeaten in [193642894.51, 226760202.29]
    dalek_G2::fixed_base_exp_window_table.push_back(193642895);
    // window 22 is unbeaten in [226760202.29, inf]
    dalek_G2::fixed_base_exp_window_table.push_back(226760202);

    /* pairing parameters */

    dalek_ate_loop_count = bigint_q("4492509698523932320491110403");
    dalek_final_exponent = bigint<6*dalek_q_limbs>("36943107177961694649618797346446870138748651578611748415128207429491593976636391130175425245705674550269561361208979548749447898941828686017765730419416875539615941651269793928962468899856083169227457503942470721108165443528513330156264699608120624990672333642644221591552000");
    dalek_final_exponent_last_chunk_abs_of_w0 = bigint_q("17970038794095729281964441603");
    dalek_final_exponent_last_chunk_is_w0_neg = true;
    dalek_final_exponent_last_chunk_w1 = bigint_q("4");

}
} // libff
