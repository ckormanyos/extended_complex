///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2016 - 2025.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#if !defined(BOOST_MATH_STANDALONE)
#define BOOST_MATH_STANDALONE
#endif

#if !defined(BOOST_MP_STANDALONE)
#define BOOST_MP_STANDALONE
#endif

#include <extended_complex.h>
#include <stopwatch.h>
#include <util.h>

#if defined(EXTENDED_COMPLEX_USE_CPP_BIN_FLOAT)
#include <boost/multiprecision/cpp_bin_float.hpp>
#else
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif
#include <boost/core/lightweight_test.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// cd /mnt/c/Users/ckorm/Documents/Ks/PC_Software/NumericalPrograms/ExtendedNumberTypes/extended_complex
// g++ -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow -O3 -std=c++23 -I. -I/mnt/c/ChrisGitRepos/modular_boost/multiprecision/include -I/mnt/c/ChrisGitRepos/modular_boost/math/include -I/mnt/c/boost/boost_1_88_0 example/example023_riemann_zeta_z.cpp example/example023a_riemann_zeta_zeros.cpp test.cpp -o test

// cd /mnt/c/Users/ckorm/Documents/Ks/PC_Software/NumericalPrograms/ExtendedNumberTypes/extended_complex/.gcov/make
// make prepare -f make_gcov_01_generic.gmk MY_ALL_COV=0 MY_BOOST_ROOT=/mnt/c/boost/boost_1_88_0 MY_CC=g++
// make gcov -f make_gcov_01_generic.gmk --jobs=2 MY_ALL_COV=0 MY_BOOST_ROOT=/mnt/c/boost/boost_1_88_0 MY_CC=g++

namespace local
{
  template<typename FloatType>
  auto my_lexical_cast(const char* p_str) -> FloatType
  {
    std::stringstream strm;

    strm << p_str;

    using float_type = FloatType;

    float_type flt;

    strm >> flt;

    return flt;
  };

  template<typename float_type>
  auto test() -> bool
  {
    using std::acos;
    using std::acosh;
    using std::asin;
    using std::atan;
    using std::atanh;
    using std::cos;
    using std::cosh;
    using std::exp;
    using std::log;
    using std::log10;
    using std::pow;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    using std::tan;
    using std::tanh;

    typedef typename extended_complex::complex<float_type> complex_type;

    const auto str_tol =
      std::string
      {
          "0."
        + std::string(static_cast<std::size_t>(std::numeric_limits<float_type>::digits10 - 2), '0')
        + std::string(static_cast<std::size_t>(UINT8_C(2)), '9')
      };

    const float_type tol { my_lexical_cast<float_type>(str_tol.c_str()) * 10 };

    std::cout << "Testing with tolerance: " << tol << std::endl;

    const complex_type val_z1(float_type(12U) / 10U, float_type(34U) / 10U);
    const complex_type val_z2(float_type(56U) / 10U, float_type(78U) / 10U);
    const complex_type val_im(float_type(0U)       , float_type(1U));

    // See also, for example, numerical evaluation at Wolfram's Alpha.

    std::vector<complex_type> results { };

    results.reserve(std::size_t { UINT8_C(64) });

    results.push_back(val_z1 / val_z2);                                            // N[((12/10) + ((34 I)/10)) / ((56/10) + ((78 I)/10)), 100]
    results.push_back(complex_type(val_z1)); results.back() /= val_z2;             // Same as above.
    results.push_back(val_z1 / (val_im * val_z2));                                 // N[((12/10) + ((34 I)/10)) / ((-78/10) + ((56 I)/10)), 100]
    results.push_back(complex_type(val_z1)); results.back() /= (val_im * val_z2);  // Same as above.
    results.push_back(val_z1.real() / val_z2);                                     // N[((12/10) / ((56/10) + ((78 I)/10)), 100]
    results.push_back(val_z1.real() / (val_im * val_z2));                          // N[((12/10) / ((-78/10) + ((56 I)/10)), 100]
    results.push_back(sqrt ( val_z1));                                             // N[Sqrt[(12/10) + ((34 I)/10)], 100]
    results.push_back(sqrt (-val_z1));                                             // N[Sqrt[(-12/10) - ((34 I)/10)], 100]
    results.push_back(sin  ( val_z1));                                             // N[Sin[(12/10) + ((34 I)/10)], 100]
    results.push_back(sinh ( val_z1));                                             // N[Sinh[(12/10) + ((34 I)/10)], 100]
    results.push_back(cosh ( val_z1));                                             // N[Cosh[(12/10) + ((34 I)/10)], 100]
    results.push_back(log  ( val_z1));                                             // N[Log[(12/10) + ((34 I)/10)], 100]
    results.push_back(asin ( val_z1));                                             // N[ArcSin[(12/10) + ((34 I)/10)], 100]
    results.push_back(acos ( val_z1));                                             // N[ArcCos[(12/10) + ((34 I)/10)], 100]
    results.push_back(atan ( val_z1));                                             // N[ArcTan[(12/10) + ((34 I)/10)], 100]
    results.push_back(acosh( val_z1));                                             // N[ArcCosh[(12/10) + ((34 I)/10)], 100]
    results.push_back(atanh( val_z1));                                             // N[ArcTanh[(12/10) + ((34 I)/10)], 100]
    results.push_back(exp  ( val_z1));                                             // N[Exp[(12/10) + ((34 I)/10)], 100]
    results.push_back(pow  ( val_z1, 17));                                         // N[((12/10) + ((34 I)/10)) ^ 17, 100]
    results.push_back(pow  ( val_z1, val_z2));                                     // N[((12/10) + ((34 I)/10)) ^ ((56/10) + ((78 I)/10)), 100]
    results.push_back(pow  ( val_z1.real(), val_z2));                              // N[(12/10)^((56/10) + ((78 I)/10)), 100]
    results.push_back(pow  ( val_z1, -17));                                        // N[((12/10) + ((34 I)/10)) ^ -17, 100]
    results.push_back(cos  ( val_z1));                                             // N[Cos[(12/10) + ((34 I)/10)], 100]
    results.push_back(pow  ( val_z1, 0));                                          // 1
    results.push_back(pow  ( val_z1, 1));                                          // 12/10 + (34 I)/10
    results.push_back(pow  ( val_z1, 2));                                          // N[((12/10) + ((34 I)/10)) ^ 2, 100]
    results.push_back(pow  ( val_z1, 3));                                          // N[((12/10) + ((34 I)/10)) ^ 3, 100]
    results.push_back(pow  ( val_z1, 4));                                          // N[((12/10) + ((34 I)/10)) ^ 4, 100]
    results.push_back(log10( val_z1));                                             // N[Log[10, (12/10) + ((34 I)/10)], 100]
    results.push_back(cosh ( val_z1));                                             // N[Cosh[(12/10) + ((34 I)/10)], 100]
    results.push_back(sinh ( val_z1));                                             // N[Sinh[(12/10) + ((34 I)/10)], 100]
    results.push_back(tanh ( val_z1));                                             // N[Tanh[(12/10) + ((34 I)/10)], 100]
    results.push_back(tan  ( val_z1));                                             // N[Tan[(12/10) + ((34 I)/10)], 100]
    using extended_complex::polar;                                                 // N[(12/10) Cos[34/10], 100]
    results.push_back(polar(real(val_z1), imag(val_z1)));                          // N[(12/10) Sin[34/10], 100]
    results.push_back(asinh( val_z1));                                             // N[ArcSinh[(12/10) + ((34 I)/10)], 100]
    results.push_back(pow  ( val_z1, real(val_z2)));                               // N[((12/10) + ((34  I)/10))^(56/10), 100]

    bool result_is_ok { true };

    // Ensure that array-oriented access works.
    {
      const bool result_access_00_is_ok = ((reinterpret_cast<const float_type*>(&val_z1))[size_t { UINT8_C(0) }] == float_type(12U) / 10U);
      const bool result_access_01_is_ok = ((reinterpret_cast<const float_type*>(&val_z1))[size_t { UINT8_C(1) }] == float_type(34U) / 10U);

      result_is_ok = (result_access_00_is_ok && result_is_ok);
      result_is_ok = (result_access_01_is_ok && result_is_ok);
    }

    // Ensure that I/O-streaming works.
    {
      std::stringstream strm;

      strm << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10))
           << results[6U];

      complex_type ctrl;

      strm >> ctrl;

      const auto result_strm_is_ok = (   util::is_close_fraction(results[6U].real(), ctrl.real(), tol)
                                      && util::is_close_fraction(results[6U].imag(), ctrl.imag(), tol));

      result_is_ok = (result_strm_is_ok && result_is_ok);
    }

    // Ensure that various forms of I/O-streaming work.
    {
      {
        std::stringstream strm;

        strm << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10))
             << "(1.25)";

        complex_type ctrl { };

        strm >> ctrl;

        using std::fpclassify;

        const auto result_strm_is_ok = ((ctrl.real() == float_type { 1.25L }) && (fpclassify(ctrl.imag()) == FP_ZERO));

        result_is_ok = (result_strm_is_ok && result_is_ok);
      }

      {
        std::stringstream strm;

        strm << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10))
             << "1.25";

        complex_type ctrl { };

        strm >> ctrl;

        using std::fpclassify;

        const auto result_strm_is_ok = ((ctrl.real() == float_type { 1.25L }) && (fpclassify(ctrl.imag()) == FP_ZERO));

        result_is_ok = (result_strm_is_ok && result_is_ok);
      }

      {
        std::stringstream strm;

        strm << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10))
             << "(x1.25)";

        complex_type ctrl { };

        strm >> ctrl;

        const auto result_strm_has_fail_is_ok = (strm.rdstate() == std::ios::failbit);;

        result_is_ok = (result_strm_has_fail_is_ok && result_is_ok);
      }
    }

    // Ensure that division with self is OK.
    {
      for(float_type flt = float_type { 1.25L }; flt < float_type { 3.0L }; flt += float_type { 0.25L })
      {
        complex_type cpx { flt, float_type { flt * float_type { 1.25L } } };

        #if defined(__clang__)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wself-assign-overloaded"
        #endif

        cpx /= cpx;

        #if defined(__clang__)
        #pragma GCC diagnostic pop
        #endif

        using std::fpclassify;

        const auto result_self_div_is_ok = ((cpx.real() == float_type { 1 }) && (fpclassify(cpx.imag()) == FP_ZERO));

        result_is_ok = (result_self_div_is_ok && result_is_ok);
      }
    }

    // Ensure that division by self is OK.
    {
      for(float_type flt = float_type { 1.25L }; flt < float_type { 3.0L }; flt += float_type { 0.25L })
      {
        complex_type cpx { flt, float_type { flt * float_type { 1.25L } } };

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

    std::vector<complex_type> controls { };

    results.reserve(std::size_t { UINT8_C(64) });

    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.3605206073752711496746203904555314533622559652928416485900216919739696312364425162689804772234273319"),    my_lexical_cast<float_type>("+0.1049891540130151843817787418655097613882863340563991323210412147505422993492407809110629067245119306")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.3605206073752711496746203904555314533622559652928416485900216919739696312364425162689804772234273319"),    my_lexical_cast<float_type>("+0.1049891540130151843817787418655097613882863340563991323210412147505422993492407809110629067245119306")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.1049891540130151843817787418655097613882863340563991323210412147505422993492407809110629067245119306"),    my_lexical_cast<float_type>("-0.3605206073752711496746203904555314533622559652928416485900216919739696312364425162689804772234273319")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.1049891540130151843817787418655097613882863340563991323210412147505422993492407809110629067245119306"),    my_lexical_cast<float_type>("-0.3605206073752711496746203904555314533622559652928416485900216919739696312364425162689804772234273319")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.07288503253796095444685466377440347071583514099783080260303687635574837310195227765726681127982646421"),   my_lexical_cast<float_type>("-0.10151843817787418655097613882863340563991323210412147505422993492407809110629067245119305856832971800")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-0.10151843817787418655097613882863340563991323210412147505422993492407809110629067245119305856832971800"),   my_lexical_cast<float_type>("-0.07288503253796095444685466377440347071583514099783080260303687635574837310195227765726681127982646421")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1.5500889128472581416161256546038815669761567486848749301860666965618993040312647033986371788677357208"),    my_lexical_cast<float_type>("+1.096711282759503047577277387056220643003106823143745046422869808875853261131777962620301480493467395")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1.096711282759503047577277387056220643003106823143745046422869808875853261131777962620301480493467395"),     my_lexical_cast<float_type>("-1.550088912847258141616125654603881566976156748684874930186066696561899304031264703398637178867735721")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+13.97940880601799793712580492576613541257396172944193599059708688128463118206190268215536541838594224"),      my_lexical_cast<float_type>("+5.42281547246340124509840716106599160358961329374827042575715571177243361237429170135167564889390308")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-1.459344510181031985739679928789446132188487461323488604725673812272622166868694452733557505015403343"),     my_lexical_cast<float_type>("-0.462696919065088203665190427736980818788403809123239459086853242811288735966522197819049006036217659")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-1.750538529873144139045226521462954860931070703406867443705575327698120127693949444003491179539803540"),     my_lexical_cast<float_type>("-0.385729418228941114585783287542904778761684113049496885765699906882071623161614865310950661528173196")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1.282474678730768368026743720782659302402633972380103558209522755331732333662205089699787331720244744"),     my_lexical_cast<float_type>("+1.231503712340851938048420309342408065643217837171236736591653326549432606404929552637127722999523972")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.327743052014525194927829972510958755346463574500092232271394201982853487105798907461836716106793827"),     my_lexical_cast<float_type>("+1.990465064891068704855135027843677587369707826516050430927052768488360486375325411568355926052994795")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1.243053274780371424303491719128792686752121125187460678216078094171054716037305591852180696564264707"),     my_lexical_cast<float_type>("-1.990465064891068704855135027843677587369707826516050430927052768488360486375325411568355926052994795")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1.472098546869956240046296809042356295374792147793626859728627826033391218153982606447668498652414512"),     my_lexical_cast<float_type>("+0.265217990171315665670057272294610940896446740868740866767584213220467425714221740314551620202361000")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1.990465064891068704855135027843677587369707826516050430927052768488360486375325411568355926052994795"),     my_lexical_cast<float_type>("+1.243053274780371424303491719128792686752121125187460678216078094171054716037305591852180696564264707")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.0865690591794584441708728351688739957204743888691393886743981396542994588169512931672375066385325971"),    my_lexical_cast<float_type>("+1.3130218230654070842957152910299990808349983340931830855923498117237087784015805687823303378344863815")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-3.209883040054176124784906450252400993119558164730356048431249139970742294562643896737048684555206883"),     my_lexical_cast<float_type>("-0.848426337294029318250973715279885597550087922172736344852553149693360359128137063129999667564390855")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-1449162978.643610460951674880000000000000000000000000000000000000000000000000000000000000000000000000"),     my_lexical_cast<float_type>("+2559363706.218816517828771840000000000000000000000000000000000000000000000000000000000000000000000000")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-0.03277613870122620601990858385164868868755535963372013573012556184586155243067581560513902047571175876"),   my_lexical_cast<float_type>("-0.08229096285844296094766104540456274850393339107196281307901532754610012233461478959682645571793968423")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.411234943477464115466545795217592784968613499972731227382657362718492707513495252941835813107553442"),     my_lexical_cast<float_type>("+2.745341999926603737618437066482640101524732307796305942046035072295581269096378050886721641340275877")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-1.675252371850406899814405010415989537774631374868045423015789145275321916111080349152932009587122213E-10"), my_lexical_cast<float_type>("-2.958659710782851665531275356243445397830042290616607863429795098837043056911533376739287222040650405E-10")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+5.43490853562576882835675504196105916449698460487268936279753180054849893360752143763634962250861355"),      my_lexical_cast<float_type>("-13.94830361398843812562310928749334812612948167559995042721295554487075295893311479161719058895913731")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1"),                                                                                                         my_lexical_cast<float_type>("+0")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1.2"),                                                                                                       my_lexical_cast<float_type>("+3.4")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-10.12"),                                                                                                     my_lexical_cast<float_type>("+8.16")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-39.888"),                                                                                                    my_lexical_cast<float_type>("-24.616")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+35.8288"),                                                                                                   my_lexical_cast<float_type>("-165.1584")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.5569716761534183846032525789711642154148645941935341359005954874987765458150971204038237271294498298"),    my_lexical_cast<float_type>("+0.5348352667130015664636074917527317522518834314413225905061936414812229669489254880191329916418075640")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-1.750538529873144139045226521462954860931070703406867443705575327698120127693949444003491179539803540"),     my_lexical_cast<float_type>("-0.385729418228941114585783287542904778761684113049496885765699906882071623161614865310950661528173196")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-1.459344510181031985739679928789446132188487461323488604725673812272622166868694452733557505015403343"),     my_lexical_cast<float_type>("-0.462696919065088203665190427736980818788403809123239459086853242811288735966522197819049006036217659")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.8505969575493737670772866494756194011367091775650155712274948188544893109556716018431217661909071889"),    my_lexical_cast<float_type>("+0.0768887100657045933256083080177585250084258780189064853743639639730289473609169028321092839605973771")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+0.0015071018758057830933874042178907505360991883301966324249948453198296050482237926843006572695221630"),    my_lexical_cast<float_type>("+1.0016427969891410443314045044734637202892888728590729858973991234850406327825813902493396819017587343")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("-1.160157831095353217138641847718832693435133312218689112813374970441475461642610761182658854558181536"),     my_lexical_cast<float_type>("-0.3066493224321975830998829152364868909731044553212128930029322191624776241668759313403878234532577257")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1.960545624274756532757863147926614306606344023611895252748744677291286163242757958621205801565261868"),     my_lexical_cast<float_type>("+1.218868916639890129907167289780557894460959308403278313426586722860313100564201336738382323265397208")));
    controls.push_back(complex_type(my_lexical_cast<float_type>("+1075.6805191639861506830414639044966711703716516323548185326564518164749257825424444319746454671929861"),    my_lexical_cast<float_type>("+757.0056142100587075756938634657815779527417605180144717715079634951737840220579027084893574878213686")));

    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 0) }].real(), controls[std::size_t { UINT8_C( 0) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 0) }].imag(), controls[std::size_t { UINT8_C( 0) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 1) }].real(), controls[std::size_t { UINT8_C( 1) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 1) }].imag(), controls[std::size_t { UINT8_C( 1) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 2) }].real(), controls[std::size_t { UINT8_C( 2) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 2) }].imag(), controls[std::size_t { UINT8_C( 2) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 3) }].real(), controls[std::size_t { UINT8_C( 3) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 3) }].imag(), controls[std::size_t { UINT8_C( 3) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 4) }].real(), controls[std::size_t { UINT8_C( 4) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 4) }].imag(), controls[std::size_t { UINT8_C( 4) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 5) }].real(), controls[std::size_t { UINT8_C( 5) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 5) }].imag(), controls[std::size_t { UINT8_C( 5) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 6) }].real(), controls[std::size_t { UINT8_C( 6) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 6) }].imag(), controls[std::size_t { UINT8_C( 6) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 7) }].real(), controls[std::size_t { UINT8_C( 7) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 7) }].imag(), controls[std::size_t { UINT8_C( 7) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 8) }].real(), controls[std::size_t { UINT8_C( 8) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 8) }].imag(), controls[std::size_t { UINT8_C( 8) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 9) }].real(), controls[std::size_t { UINT8_C( 9) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C( 9) }].imag(), controls[std::size_t { UINT8_C( 9) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(10) }].real(), controls[std::size_t { UINT8_C(10) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(10) }].imag(), controls[std::size_t { UINT8_C(10) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(11) }].real(), controls[std::size_t { UINT8_C(11) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(11) }].imag(), controls[std::size_t { UINT8_C(11) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(12) }].real(), controls[std::size_t { UINT8_C(12) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(12) }].imag(), controls[std::size_t { UINT8_C(12) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(13) }].real(), controls[std::size_t { UINT8_C(13) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(13) }].imag(), controls[std::size_t { UINT8_C(13) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(14) }].real(), controls[std::size_t { UINT8_C(14) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(14) }].imag(), controls[std::size_t { UINT8_C(14) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(15) }].real(), controls[std::size_t { UINT8_C(15) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(15) }].imag(), controls[std::size_t { UINT8_C(15) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(16) }].real(), controls[std::size_t { UINT8_C(16) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(16) }].imag(), controls[std::size_t { UINT8_C(16) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(17) }].real(), controls[std::size_t { UINT8_C(17) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(17) }].imag(), controls[std::size_t { UINT8_C(17) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(18) }].real(), controls[std::size_t { UINT8_C(18) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(18) }].imag(), controls[std::size_t { UINT8_C(18) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(19) }].real(), controls[std::size_t { UINT8_C(19) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(19) }].imag(), controls[std::size_t { UINT8_C(19) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(20) }].real(), controls[std::size_t { UINT8_C(20) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(20) }].imag(), controls[std::size_t { UINT8_C(20) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(21) }].real(), controls[std::size_t { UINT8_C(21) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(21) }].imag(), controls[std::size_t { UINT8_C(21) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(22) }].real(), controls[std::size_t { UINT8_C(22) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(22) }].imag(), controls[std::size_t { UINT8_C(22) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(23) }].real(), controls[std::size_t { UINT8_C(23) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(23) }].imag(), controls[std::size_t { UINT8_C(23) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(24) }].real(), controls[std::size_t { UINT8_C(24) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(24) }].imag(), controls[std::size_t { UINT8_C(24) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(25) }].real(), controls[std::size_t { UINT8_C(25) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(25) }].imag(), controls[std::size_t { UINT8_C(25) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(26) }].real(), controls[std::size_t { UINT8_C(26) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(26) }].imag(), controls[std::size_t { UINT8_C(26) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(27) }].real(), controls[std::size_t { UINT8_C(27) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(27) }].imag(), controls[std::size_t { UINT8_C(27) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(28) }].real(), controls[std::size_t { UINT8_C(28) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(28) }].imag(), controls[std::size_t { UINT8_C(28) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(29) }].real(), controls[std::size_t { UINT8_C(29) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(29) }].imag(), controls[std::size_t { UINT8_C(29) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(30) }].real(), controls[std::size_t { UINT8_C(30) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(30) }].imag(), controls[std::size_t { UINT8_C(30) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(31) }].real(), controls[std::size_t { UINT8_C(31) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(31) }].imag(), controls[std::size_t { UINT8_C(31) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(32) }].real(), controls[std::size_t { UINT8_C(32) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(32) }].imag(), controls[std::size_t { UINT8_C(32) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(33) }].real(), controls[std::size_t { UINT8_C(33) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(33) }].imag(), controls[std::size_t { UINT8_C(33) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(34) }].real(), controls[std::size_t { UINT8_C(34) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(34) }].imag(), controls[std::size_t { UINT8_C(34) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }
    { const auto result_real_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(35) }].real(), controls[std::size_t { UINT8_C(35) }].real(), tol); const auto result_imag_is_ok = util::is_close_fraction(results[std::size_t { UINT8_C(35) }].imag(), controls[std::size_t { UINT8_C(35) }].imag(), tol); BOOST_TEST(result_real_is_ok); BOOST_TEST(result_imag_is_ok); result_is_ok = (result_real_is_ok && result_imag_is_ok && result_is_ok); std::stringstream strm { }; strm << "result_real_is_ok: " << std::boolalpha << result_real_is_ok << ", result_imag_is_ok: " << std::boolalpha << result_imag_is_ok; std::cout << strm.str() << std::endl; }

    return result_is_ok;
  }
}

extern auto example023_riemann_zeta_z     () -> bool;
extern auto example023a_riemann_zeta_zeros() -> bool;

auto main() -> int
{
  #if defined(EXTENDED_COMPLEX_USE_CPP_BIN_FLOAT)
  using local_mp_type = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<100>, boost::multiprecision::et_off>;
  #else
  using local_mp_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>, boost::multiprecision::et_off>;
  #endif

  using stopwatch_type = concurrency::stopwatch;

  stopwatch_type my_stopwatch { };

  const auto result_flt_________is_ok = local::test<float>();
  const auto result_dbl_________is_ok = local::test<double>();
  const auto result_ldbl________is_ok = local::test<long double>();
  const auto result_mp__________is_ok = local::test<local_mp_type>();
  const auto result_example023__is_ok = ::example023_riemann_zeta_z();
  const auto result_example023a_is_ok = ::example023a_riemann_zeta_zeros();

  const auto execution_time = stopwatch_type::elapsed_time<float>(my_stopwatch);

  const bool result_is_ok
  {
       result_flt_________is_ok
    && result_dbl_________is_ok
    && result_ldbl________is_ok
    && result_mp__________is_ok
    && result_example023__is_ok
    && result_example023a_is_ok
  };

  {
    std::stringstream strm { };

    strm << "\nresult_is_ok: "
         << std::boolalpha
         << result_is_ok
         << ", time: "
         << std::fixed
         << std::setprecision(1)
         << execution_time
         << "s";

    std::cout << strm.str() << std::endl;
  }

  BOOST_TEST(result_is_ok);

  BOOST_TEST(result_flt_________is_ok);
  BOOST_TEST(result_dbl_________is_ok);
  BOOST_TEST(result_ldbl________is_ok);
  BOOST_TEST(result_mp__________is_ok);
  BOOST_TEST(result_example023__is_ok);
  BOOST_TEST(result_example023a_is_ok);

  return boost::report_errors();
}
