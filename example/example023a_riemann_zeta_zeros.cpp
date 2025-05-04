///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2024 - 2025.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <example/zeta_detail.h>
#include <extended_complex.h>

#include <boost/math/tools/toms748_solve.hpp>
#if defined(EXTENDED_COMPLEX_USE_CPP_BIN_FLOAT)
#include <boost/multiprecision/cpp_bin_float.hpp>
#else
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#if defined(__GNUC__)
#define EXTENDED_COMPLEX_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define EXTENDED_COMPLEX_NOINLINE __declspec(noinline)
#else
#define EXTENDED_COMPLEX_NOINLINE
#endif

namespace local
{
  namespace detail
  {
    #if !defined(EXTENDED_COMPLEX_REDUCE_TEST_DEPTH)
    constexpr unsigned multiprecision_digits10 { static_cast<unsigned>(UINT16_C(501)) };
    #else
    constexpr unsigned multiprecision_digits10 { static_cast<unsigned>(UINT16_C(101)) };
    #endif

    #if defined(EXTENDED_COMPLEX_USE_CPP_BIN_FLOAT)
    using multiprecision_float_type =
      boost::multiprecision::number<boost::multiprecision::cpp_bin_float<multiprecision_digits10>,
                                    boost::multiprecision::et_off>;
    #else
    using multiprecision_float_type =
      boost::multiprecision::number<boost::multiprecision::cpp_dec_float<multiprecision_digits10>,
                                    boost::multiprecision::et_off>;
    #endif
  }

  using complex_type = extended_complex::complex<detail::multiprecision_float_type>;

  using real_type = typename complex_type::value_type;

  using bracket_type = std::pair<real_type, real_type>;

  EXTENDED_COMPLEX_NOINLINE static auto find_riemann_root(const bracket_type& bt_val) -> complex_type
  {
    auto tol
    {
      [](const real_type& a, const real_type& b)
      {
        using std::fabs;

        const real_type delta { fabs(1 - (a / b)) };

        return (delta < real_type { std::numeric_limits<real_type>::epsilon() * static_cast<unsigned>(UINT8_C(64)) });
      }
    };

    static const real_type my_half { 0.5F };

    auto
      my_riemann_function
      {
        [](const real_type& y) -> real_type
        {
          return riemann_zeta(complex_type { my_half, y }).real();
        }
      };

    std::uintmax_t max_iter { static_cast<std::uintmax_t>(UINT8_C(64)) };

    const std::pair<real_type, real_type>
      result
      {
        boost::math::tools::toms748_solve(my_riemann_function, bt_val.first, bt_val.second, tol, max_iter)
      };

    {
      std::stringstream strm { };

      strm << std::setprecision(std::numeric_limits<real_type>::digits10) << result.first;

      std::cout << strm.str() << std::endl;
    }

    return { my_half, result.first };
  }
} // namespace local

auto example023a_riemann_zeta_zeros() -> bool
{
  using local::real_type;
  using local::bracket_type;

  const bracket_type
    bt_val0
    {
      real_type { static_cast<unsigned>(UINT16_C(141)) } / static_cast<unsigned>(UINT8_C(10)),
      real_type { static_cast<unsigned>(UINT16_C(142)) } / static_cast<unsigned>(UINT8_C(10)),
    };

  const bracket_type
    bt_val1
    {
      real_type { static_cast<unsigned>(UINT16_C(210)) } / static_cast<unsigned>(UINT8_C(10)),
      real_type { static_cast<unsigned>(UINT16_C(211)) } / static_cast<unsigned>(UINT8_C(10)),
    };

  const bracket_type
    bt_val2
    {
      real_type { static_cast<unsigned>(UINT16_C(250)) } / static_cast<unsigned>(UINT8_C(10)),
      real_type { static_cast<unsigned>(UINT16_C(251)) } / static_cast<unsigned>(UINT8_C(10)),
    };

  const bracket_type
    bt_val3
    {
      real_type { static_cast<unsigned>(UINT16_C(304)) } / static_cast<unsigned>(UINT8_C(10)),
      real_type { static_cast<unsigned>(UINT16_C(305)) } / static_cast<unsigned>(UINT8_C(10)),
    };

  const bracket_type
    bt_val4
    {
      real_type { static_cast<unsigned>(UINT16_C(329)) } / static_cast<unsigned>(UINT8_C(10)),
      real_type { static_cast<unsigned>(UINT16_C(330)) } / static_cast<unsigned>(UINT8_C(10)),
    };

  const bracket_type
    bt_val5
    {
      real_type { static_cast<unsigned>(UINT16_C(375)) } / static_cast<unsigned>(UINT8_C(10)),
      real_type { static_cast<unsigned>(UINT16_C(376)) } / static_cast<unsigned>(UINT8_C(10)),
    };

  const bracket_type
    bt_val6
    {
      real_type { static_cast<unsigned>(UINT16_C(409)) } / static_cast<unsigned>(UINT8_C(10)),
      real_type { static_cast<unsigned>(UINT16_C(410)) } / static_cast<unsigned>(UINT8_C(10)),
    };

  // rz0: 14.1347251417346937904572519835624702707842571156992431756855674601499634298092567649490103931715610127792029715487974367661426914698822545825053632394471377804133812372059705496219558658602005555667258360107737002054109826615075427805174425913062544819786510723049387256297383215774203952157256748093321400349904680343462673144209203773854871413783173563969953654281130796805314916885290678208229804926433866673462332007875876179200560486805435680144442465106559756866590322868651054485944432062407271
  // rz1: 21.0220396387715549926284795938969027773343405249027817546295204035875985860688907997136585141801514195337254736424758913838650686037313212621188216243757416692565447118440711940313067256462277926148873374355520591473971328226624707890767538144407264668419060771275698340545140284399232225367882682361112892700575856532731588666042140009071151080090069720027998711017584751963221649686590057481124793869163835183723427807344902391010385045756412159583999210016218346691131587217480571703157935817977250
  // rz2: 25.0108575801456887632137909925628218186595496725579966724965420067450920984416442778402382245580624407504710461490557783782998515227308011881339335826716895872251698104387355129284937271919946229759126754786966288568077350700399577231140232842768736693998732195864877522500991924534749762085766123345997354435583675313812659977645290374484969947911378977220661993071899723225497322716300515916192127977408766000672914983081279306670273508495160019846705424694917966952255141793196653912734145216731598
  // rz3: 30.4248761258595132103118975305840913201815600237154401809621460369933293893332779202905842939020891106309917115273954991176332266711863193918072259567142433411559068546813655807241734984472495931904081163231501970234848416302214009856207397183920181330218680632982257197522500237468561369747124964426229779245040574906715345727886515065160832468797062817781045777722587891923738629001127603097356808904925300646128927275309194479020035893898194274955113239173842716381084004992111980069243871887296959
  // rz4: 32.9350615877391896906623689640749034888127156035170390092800034407848156086305510059388484961353487245496025252805975815135792377828577545060376530114726821098252727136594781660791865078811703538367654746017385481206517878865964665947287871860279716580426776485440666929093931931156455083917513007902799895626566839200248741651699908688450128764042509560641448930237747970532527824099477517659446074468098740670653382909589156142048002241984083751415844352378065847539082079012460705045601001022574389
  // rz5: 37.5861781588256712572177634807053328214055973508307932183330011136221490896185372647303291049458238034745777477461922353799409650222639362856248220474836880890036678903545375537779031909056440993758387212756343126430515074908971220774678015512043279929866514323315465488505867950603413625492599668820961532811560379952622032357848783110737761462301491055793790156771634112512653712956501767448864521498327020119844943153262080331531776074647950750405608208302136339508622882205513357695558035937464812
  // rz6: 40.9187190121474951873981269146332543957261659627772795361613036672532805287200712829960037198895468755036655196769454679449171364948985991756534367533105048091521687544875227792200700073971948417805644552230105285546119797579028708944367274046299954096047330409695262215628459484299449278461845787877616326789854362772538625569030381839732705149501966655277103497437488166959905285910304820629636478328802052219973188256580123512214090671916895581038535743410889979247590603457572189960393020947459656

  using local::complex_type;
  using local::find_riemann_root;

  // TODO: Find the roots with parallel-for, but also notice that
  // riemann_zeta(...) might not yet be strictly parallelizable.

  const complex_type rz0 { find_riemann_root(bt_val0) };
  const complex_type rz1 { find_riemann_root(bt_val1) };
  const complex_type rz2 { find_riemann_root(bt_val2) };
  const complex_type rz3 { find_riemann_root(bt_val3) };
  const complex_type rz4 { find_riemann_root(bt_val4) };
  const complex_type rz5 { find_riemann_root(bt_val5) };
  const complex_type rz6 { find_riemann_root(bt_val6) };

  const real_type my_real_tol { std::numeric_limits<real_type>::epsilon() * static_cast<unsigned>(UINT32_C(0x100000)) };

  const bool result_rz0_is_ok { (abs(riemann_zeta(rz0)) < my_real_tol) };
  const bool result_rz1_is_ok { (abs(riemann_zeta(rz1)) < my_real_tol) };
  const bool result_rz2_is_ok { (abs(riemann_zeta(rz2)) < my_real_tol) };
  const bool result_rz3_is_ok { (abs(riemann_zeta(rz3)) < my_real_tol) };
  const bool result_rz4_is_ok { (abs(riemann_zeta(rz4)) < my_real_tol) };
  const bool result_rz5_is_ok { (abs(riemann_zeta(rz5)) < my_real_tol) };
  const bool result_rz6_is_ok { (abs(riemann_zeta(rz6)) < my_real_tol) };

  const bool
    result_is_ok
    {
         result_rz0_is_ok
      && result_rz1_is_ok
      && result_rz2_is_ok
      && result_rz3_is_ok
      && result_rz4_is_ok
      && result_rz5_is_ok
      && result_rz6_is_ok
    };

  return result_is_ok;
}
