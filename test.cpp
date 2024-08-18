///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2016 - 2024.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/core/lightweight_test.hpp>

#include <extended_complex.h>

// cd /mnt/c/Users/ckorm/Documents/Ks/PC_Software/NumericalPrograms/ExtendedNumberTypes/extended_complex
// g++ -Werror -Wall -Wextra -Wno-unused-parameter -O3 -std=c++20 -I. -I/mnt/c/ChrisGitRepos/modular_boost/multiprecision/include -I/mnt/c/boost/boost_1_85_0 example/example023_riemann_zeta_z.cpp example/example023a_riemann_zeta_zeros.cpp test.cpp -o test.exe

// cd /mnt/c/Users/ckorm/Documents/Ks/PC_Software/NumericalPrograms/ExtendedNumberTypes/extended_complex/.gcov/make
// make prepare -f make_gcov_01_generic.gmk MY_ALL_COV=0 MY_BOOST_ROOT=/mnt/c/boost/boost_1_85_0 MY_CC=g++
// make gcov -f make_gcov_01_generic.gmk --jobs=2 MY_ALL_COV=0 MY_BOOST_ROOT=/mnt/c/boost/boost_1_85_0 MY_CC=g++


namespace local
{
  template<typename NumericType>
  auto is_close_fraction(const NumericType& a,
                         const NumericType& b,
                         const NumericType& tol) noexcept -> bool
  {
    using std::fabs;

    auto result_is_ok = bool { };

    if(b == static_cast<NumericType>(0))
    {
      result_is_ok = (fabs(a - b) < tol);
    }
    else
    {
      const auto delta = fabs(1 - (a / b));

      result_is_ok = (delta < tol);
    }

    return result_is_ok;
  }

  template<typename FloatType>
  auto my_lexical_cast(const char* p_str) -> FloatType
  {
    std::stringstream strm;

    strm << p_str;

    using float_type = FloatType;

    float_type flt;

    strm >> flt;

    return std::move(static_cast<float_type&&>(flt));
  };

  template<typename float_type>
  auto test() -> bool
  {
    using std::acos;
    using std::asin;
    using std::atan;
    using std::atanh;
    using std::cos;
    using std::exp;
    using std::log;
    using std::log10;
    using std::pow;
    using std::sin;
    using std::sinh;
    using std::sqrt;

    typedef typename extended_complex::complex<float_type> complex_type;

    const auto str_tol =
      std::string
      {
          "0."
        + std::string(static_cast<std::size_t>(std::numeric_limits<float_type>::digits10 - 2), '0')
        + std::string(static_cast<std::size_t>(UINT8_C(2)), '9')
      };

    const auto tol = my_lexical_cast<float_type>(str_tol.c_str());

    std::cout << "Testing with tolerance: " << tol << std::endl;

    const complex_type val_z1(float_type(12U) / 10U, float_type(34U) / 10U);
    const complex_type val_z2(float_type(56U) / 10U, float_type(78U) / 10U);
    const complex_type val_im(float_type(0U)       , float_type(1U));

    // See also, for example, numerical evaluation at Wolfram's Alpha.

    const complex_type result_01 = val_z1 / val_z2;                                       // N[((12/10) + ((34 I)/10)) / ((56/10) + ((78 I)/10)), 100]
          complex_type result_02 = complex_type(val_z1); result_02 /= val_z2;             // Same as above.
    const complex_type result_03 = val_z1 / (val_im * val_z2);                            // N[((12/10) + ((34 I)/10)) / ((-78/10) + ((56 I)/10)), 100]
          complex_type result_04 = complex_type(val_z1); result_04 /= (val_im * val_z2);  // Same as above.
    const complex_type result_05 = val_z1.real() / val_z2;                                // N[((12/10) / ((56/10) + ((78 I)/10)), 100]
    const complex_type result_06 = val_z1.real() / (val_im * val_z2);                     // N[((12/10) / ((-78/10) + ((56 I)/10)), 100]
    const complex_type result_07 = sqrt ( val_z1);                                        // N[Sqrt[(12/10) + ((34 I)/10)], 100]
    const complex_type result_08 = sqrt (-val_z1);                                        // N[Sqrt[(-12/10) - ((34 I)/10)], 100]
    const complex_type result_09 = sin  ( val_z1);                                        // N[Sin[(12/10) + ((34 I)/10)], 100]
    const complex_type result_10 = sinh ( val_z1);                                        // N[Sinh[(12/10) + ((34 I)/10)], 100]
    const complex_type result_11 = cosh ( val_z1);                                        // N[Cosh[(12/10) + ((34 I)/10)], 100]
    const complex_type result_12 = log  ( val_z1);                                        // N[Log[(12/10) + ((34 I)/10)], 100]
    const complex_type result_13 = asin ( val_z1);                                        // N[ArcSin[(12/10) + ((34 I)/10)], 100]
    const complex_type result_14 = acos ( val_z1);                                        // N[ArcCos[(12/10) + ((34 I)/10)], 100]
    const complex_type result_15 = atan ( val_z1);                                        // N[ArcTan[(12/10) + ((34 I)/10)], 100]
    const complex_type result_16 = acosh( val_z1);                                        // N[ArcCosh[(12/10) + ((34 I)/10)], 100]
    const complex_type result_17 = atanh( val_z1);                                        // N[ArcTanh[(12/10) + ((34 I)/10)], 100]
    const complex_type result_18 = exp  ( val_z1);                                        // N[Exp[(12/10) + ((34 I)/10)], 100]
    const complex_type result_19 = pow  ( val_z1, 17);                                    // N[((12/10) + ((34 I)/10)) ^ 17, 100]
    const complex_type result_20 = pow  ( val_z1, val_z2);                                // N[((12/10) + ((34 I)/10)) ^ ((56/10) + ((78 I)/10)), 100]
    const complex_type result_21 = pow  ( val_z1.real(), val_z2);                         // N[(12/10)^((56/10) + ((78 I)/10)), 100]
    const complex_type result_22 = pow  ( val_z1, -17);                                   // N[((12/10) + ((34 I)/10)) ^ -17, 100]
    const complex_type result_23 = cos  ( val_z1);                                        // N[Cos[(12/10) + ((34 I)/10)], 100]
    const complex_type result_24 = pow  ( val_z1, 0);                                     // 1
    const complex_type result_25 = pow  ( val_z1, 1);                                     // 12/10 + (34 I)/10
    const complex_type result_26 = pow  ( val_z1, 2);                                     // N[((12/10) + ((34 I)/10)) ^ 2, 100]
    const complex_type result_27 = pow  ( val_z1, 3);                                     // N[((12/10) + ((34 I)/10)) ^ 3, 100]
    const complex_type result_28 = pow  ( val_z1, 4);                                     // N[((12/10) + ((34 I)/10)) ^ 4, 100]
    const complex_type result_29 = log10( val_z1);                                        // N[Log[10, (12/10) + ((34 I)/10)], 100]

    auto result_is_ok = true;

    // Ensure that I/O-streaming works.
    {
      std::stringstream strm;

      strm << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10))
           << result_07;

      complex_type ctrl;

      strm >> ctrl;

      const auto result_strm_is_ok = (   is_close_fraction(result_07.real(), ctrl.real(), tol)
                                      && is_close_fraction(result_07.imag(), ctrl.imag(), tol));

      result_is_ok = (result_strm_is_ok && result_is_ok);
    }

    // Ensure that all forms of I/O-streaming work.
    {
      {
        std::stringstream strm;

        strm << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10))
             << "(1.25)";

        complex_type ctrl;

        strm >> ctrl;

        using std::fpclassify;

        const auto result_strm_is_ok = (ctrl.real() == float_type { 1.25L } && fpclassify(ctrl.imag()) == FP_ZERO);

        result_is_ok = (result_strm_is_ok && result_is_ok);
      }

      {
        std::stringstream strm;

        strm << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10))
             << "1.25";

        complex_type ctrl;

        strm >> ctrl;

        using std::fpclassify;

        const auto result_strm_is_ok = (ctrl.real() == float_type { 1.25L } && fpclassify(ctrl.imag()) == FP_ZERO);

        result_is_ok = (result_strm_is_ok && result_is_ok);
      }

      {
        std::stringstream strm;

        strm << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10))
             << "(x1.25)";

        complex_type ctrl;

        strm >> ctrl;

        const auto result_strm_has_fail_is_ok = (strm.rdstate() == std::ios_base::failbit);;

        result_is_ok = (result_strm_has_fail_is_ok && result_is_ok);
      }
    }

    // Ensure that division with self is OK.
    {
      for(float_type flt = float_type { 1.25L }; flt < float_type { 3.0L }; flt += float_type { 0.25L })
      {
        complex_type cpx;

        #if defined(__clang__)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wself-assign-overloaded"
        #endif

        cpx /= cpx;

        #if defined(__clang__)
        #pragma GCC diagnostic pop
        #endif

        using std::fpclassify;

        const auto result_self_div_is_ok = ((cpx.real() == 1) && (fpclassify(cpx.imag()) == FP_ZERO));

        result_is_ok = (result_self_div_is_ok && result_is_ok);
      }
    }

    // Ensure that division by self is OK.
    {
      for(float_type flt = float_type { 1.25L }; flt < float_type { 3.0L }; flt += float_type { 0.25L })
      {
        complex_type cpx;

        #if defined(__clang__)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wself-assign-overloaded"
        #endif

        cpx -= cpx;

        #if defined(__clang__)
        #pragma GCC diagnostic pop
        #endif

        using std::fpclassify;

        const auto result_self_sub_is_ok = (((fpclassify(cpx.real()) == FP_ZERO)) && (fpclassify(cpx.imag()) == FP_ZERO));

        result_is_ok = (result_self_sub_is_ok && result_is_ok);
      }
    }

    const complex_type control_01(my_lexical_cast<float_type>(  "+0.3605206073752711496746203904555314533622559652928416485900216919739696312364425162689804772234273319"),    my_lexical_cast<float_type>( "+0.1049891540130151843817787418655097613882863340563991323210412147505422993492407809110629067245119306"));
    const complex_type control_02(my_lexical_cast<float_type>(  "+0.3605206073752711496746203904555314533622559652928416485900216919739696312364425162689804772234273319"),    my_lexical_cast<float_type>( "+0.1049891540130151843817787418655097613882863340563991323210412147505422993492407809110629067245119306"));
    const complex_type control_03(my_lexical_cast<float_type>(  "+0.1049891540130151843817787418655097613882863340563991323210412147505422993492407809110629067245119306"),    my_lexical_cast<float_type>( "-0.3605206073752711496746203904555314533622559652928416485900216919739696312364425162689804772234273319"));
    const complex_type control_04(my_lexical_cast<float_type>(  "+0.1049891540130151843817787418655097613882863340563991323210412147505422993492407809110629067245119306"),    my_lexical_cast<float_type>( "-0.3605206073752711496746203904555314533622559652928416485900216919739696312364425162689804772234273319"));
    const complex_type control_05(my_lexical_cast<float_type>(  "+0.07288503253796095444685466377440347071583514099783080260303687635574837310195227765726681127982646421"),   my_lexical_cast<float_type>( "-0.10151843817787418655097613882863340563991323210412147505422993492407809110629067245119305856832971800"));
    const complex_type control_06(my_lexical_cast<float_type>(  "-0.10151843817787418655097613882863340563991323210412147505422993492407809110629067245119305856832971800"),   my_lexical_cast<float_type>( "-0.07288503253796095444685466377440347071583514099783080260303687635574837310195227765726681127982646421"));
    const complex_type control_07(my_lexical_cast<float_type>(  "+1.5500889128472581416161256546038815669761567486848749301860666965618993040312647033986371788677357208"),    my_lexical_cast<float_type>( "+1.096711282759503047577277387056220643003106823143745046422869808875853261131777962620301480493467395"));
    const complex_type control_08(my_lexical_cast<float_type>(  "+1.096711282759503047577277387056220643003106823143745046422869808875853261131777962620301480493467395"),     my_lexical_cast<float_type>( "-1.550088912847258141616125654603881566976156748684874930186066696561899304031264703398637178867735721"));
    const complex_type control_09(my_lexical_cast<float_type>( "+13.97940880601799793712580492576613541257396172944193599059708688128463118206190268215536541838594224"),      my_lexical_cast<float_type>( "+5.42281547246340124509840716106599160358961329374827042575715571177243361237429170135167564889390308"));
    const complex_type control_10(my_lexical_cast<float_type>(  "-1.459344510181031985739679928789446132188487461323488604725673812272622166868694452733557505015403343"),     my_lexical_cast<float_type>( "-0.462696919065088203665190427736980818788403809123239459086853242811288735966522197819049006036217659"));
    const complex_type control_11(my_lexical_cast<float_type>(  "-1.750538529873144139045226521462954860931070703406867443705575327698120127693949444003491179539803540"),     my_lexical_cast<float_type>( "-0.385729418228941114585783287542904778761684113049496885765699906882071623161614865310950661528173196"));
    const complex_type control_12(my_lexical_cast<float_type>(  "+1.282474678730768368026743720782659302402633972380103558209522755331732333662205089699787331720244744"),     my_lexical_cast<float_type>( "+1.231503712340851938048420309342408065643217837171236736591653326549432606404929552637127722999523972"));
    const complex_type control_13(my_lexical_cast<float_type>(  "+0.327743052014525194927829972510958755346463574500092232271394201982853487105798907461836716106793827"),     my_lexical_cast<float_type>( "+1.990465064891068704855135027843677587369707826516050430927052768488360486375325411568355926052994795"));
    const complex_type control_14(my_lexical_cast<float_type>(  "+1.243053274780371424303491719128792686752121125187460678216078094171054716037305591852180696564264707"),     my_lexical_cast<float_type>( "-1.990465064891068704855135027843677587369707826516050430927052768488360486375325411568355926052994795"));
    const complex_type control_15(my_lexical_cast<float_type>(  "+1.472098546869956240046296809042356295374792147793626859728627826033391218153982606447668498652414512"),     my_lexical_cast<float_type>( "+0.265217990171315665670057272294610940896446740868740866767584213220467425714221740314551620202361000"));
    const complex_type control_16(my_lexical_cast<float_type>(  "+1.990465064891068704855135027843677587369707826516050430927052768488360486375325411568355926052994795"),     my_lexical_cast<float_type>( "+1.243053274780371424303491719128792686752121125187460678216078094171054716037305591852180696564264707"));
    const complex_type control_17(my_lexical_cast<float_type>(  "+0.0865690591794584441708728351688739957204743888691393886743981396542994588169512931672375066385325971"),    my_lexical_cast<float_type>( "+1.3130218230654070842957152910299990808349983340931830855923498117237087784015805687823303378344863815"));
    const complex_type control_18(my_lexical_cast<float_type>(  "-3.209883040054176124784906450252400993119558164730356048431249139970742294562643896737048684555206883"),     my_lexical_cast<float_type>( "-0.848426337294029318250973715279885597550087922172736344852553149693360359128137063129999667564390855"));
    const complex_type control_19(my_lexical_cast<float_type>(  "-1449162978.643610460951674880000000000000000000000000000000000000000000000000000000000000000000000000"),     my_lexical_cast<float_type>( "+2559363706.218816517828771840000000000000000000000000000000000000000000000000000000000000000000000000"));
    const complex_type control_20(my_lexical_cast<float_type>(  "-0.03277613870122620601990858385164868868755535963372013573012556184586155243067581560513902047571175876"),   my_lexical_cast<float_type>( "-0.08229096285844296094766104540456274850393339107196281307901532754610012233461478959682645571793968423"));
    const complex_type control_21(my_lexical_cast<float_type>(  "+0.411234943477464115466545795217592784968613499972731227382657362718492707513495252941835813107553442"),     my_lexical_cast<float_type>( "+2.745341999926603737618437066482640101524732307796305942046035072295581269096378050886721641340275877"));
    const complex_type control_22(my_lexical_cast<float_type>(  "-1.675252371850406899814405010415989537774631374868045423015789145275321916111080349152932009587122213E-10"), my_lexical_cast<float_type>( "-2.958659710782851665531275356243445397830042290616607863429795098837043056911533376739287222040650405E-10"));
    const complex_type control_23(my_lexical_cast<float_type>(  "+5.43490853562576882835675504196105916449698460487268936279753180054849893360752143763634962250861355"),      my_lexical_cast<float_type>( "-13.94830361398843812562310928749334812612948167559995042721295554487075295893311479161719058895913731"));
    const complex_type control_24(my_lexical_cast<float_type>(  "+1"),                                                                                                         my_lexical_cast<float_type>( "+0"));
    const complex_type control_25(my_lexical_cast<float_type>(  "+1.2"),                                                                                                       my_lexical_cast<float_type>( "+3.4"));
    const complex_type control_26(my_lexical_cast<float_type>(  "-10.12"),                                                                                                     my_lexical_cast<float_type>( "+8.16"));
    const complex_type control_27(my_lexical_cast<float_type>(  "-39.888"),                                                                                                    my_lexical_cast<float_type>( "-24.616"));
    const complex_type control_28(my_lexical_cast<float_type>(  "+35.8288"),                                                                                                   my_lexical_cast<float_type>( "-165.1584"));
    const complex_type control_29(my_lexical_cast<float_type>(  "+0.5569716761534183846032525789711642154148645941935341359005954874987765458150971204038237271294498298"),    my_lexical_cast<float_type>( "+0.5348352667130015664636074917527317522518834314413225905061936414812229669489254880191329916418075640"));

    { const auto result_real_is_ok = is_close_fraction(result_01.real(), control_01.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_01.imag(), control_01.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_02.real(), control_02.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_02.imag(), control_02.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_03.real(), control_03.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_03.imag(), control_03.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_04.real(), control_04.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_04.imag(), control_04.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_05.real(), control_05.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_05.imag(), control_05.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_06.real(), control_06.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_06.imag(), control_06.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_07.real(), control_07.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_07.imag(), control_07.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_08.real(), control_08.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_08.imag(), control_08.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_09.real(), control_09.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_09.imag(), control_09.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_10.real(), control_10.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_10.imag(), control_10.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_11.real(), control_11.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_11.imag(), control_11.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_12.real(), control_12.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_12.imag(), control_12.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_13.real(), control_13.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_13.imag(), control_13.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_14.real(), control_14.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_14.imag(), control_14.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_15.real(), control_15.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_15.imag(), control_15.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_16.real(), control_16.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_16.imag(), control_16.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_17.real(), control_17.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_17.imag(), control_17.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_18.real(), control_18.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_18.imag(), control_18.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_19.real(), control_19.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_19.imag(), control_19.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_20.real(), control_20.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_20.imag(), control_20.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_21.real(), control_21.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_21.imag(), control_21.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_22.real(), control_22.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_22.imag(), control_22.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_23.real(), control_23.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_23.imag(), control_23.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_24.real(), control_24.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_24.imag(), control_24.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_25.real(), control_25.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_25.imag(), control_25.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_26.real(), control_26.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_26.imag(), control_26.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_27.real(), control_27.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_27.imag(), control_27.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_28.real(), control_28.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_28.imag(), control_28.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = is_close_fraction(result_29.real(), control_29.real(), tol); const auto result_imag_is_ok = is_close_fraction(result_29.imag(), control_29.imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }

    return result_is_ok;
  }
}

extern auto example023_riemann_zeta_z     () -> bool;
extern auto example023a_riemann_zeta_zeros() -> bool;

auto main() -> int
{
  using local_mp_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>, boost::multiprecision::et_off>;

  const auto result_flt_________is_ok = local::test<float>();
  const auto result_dbl_________is_ok = local::test<double>();
  const auto result_ldbl________is_ok = local::test<long double>();
  const auto result_mp__________is_ok = local::test<local_mp_type>();
  const auto result_example023__is_ok = ::example023_riemann_zeta_z();
  const auto result_example023a_is_ok = ::example023a_riemann_zeta_zeros();

  BOOST_TEST(result_flt_________is_ok);
  BOOST_TEST(result_dbl_________is_ok);
  BOOST_TEST(result_ldbl________is_ok);
  BOOST_TEST(result_mp__________is_ok);
  BOOST_TEST(result_example023__is_ok);
  BOOST_TEST(result_example023a_is_ok);

  auto result_is_ok =
  (
       result_flt_________is_ok
    && result_dbl_________is_ok
    && result_ldbl________is_ok
    && result_mp__________is_ok
    && result_example023__is_ok
    && result_example023a_is_ok
  );

  BOOST_TEST(result_is_ok);

  return boost::report_errors();
}
