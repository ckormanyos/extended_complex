///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2024.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <example/zeta_detail.h>
#include <extended_complex.h>

#include <boost/math/tools/toms748_solve.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

namespace local
{
  using complex_type = extended_complex::complex<boost::multiprecision::number<boost::multiprecision::cpp_dec_float<static_cast<unsigned>(UINT8_C(101))>, boost::multiprecision::et_off>>;
  using real_type    = typename complex_type::value_type;

  real_type my_riemann_function(const real_type& y)
  {
    static const real_type my_half { real_type { static_cast<int>(INT8_C(1)) } / static_cast<int>(INT8_C(2)) };

    return riemann_zeta(complex_type { my_half, y }).real();
  }

  using bracket_type = std::pair<real_type, real_type>;

  auto find_riemann_root(const bracket_type& bt_val) -> complex_type
  {
    auto tol =
      [](const real_type& a, const real_type& b)
      {
        const real_type delta { fabs(1 - (a / b)) };

        static const real_type my_tol { std::numeric_limits<real_type>::epsilon() * 64 };

        return (delta < my_tol);
      };

    std::uintmax_t max_iter { static_cast<std::uintmax_t>(UINT8_C(32)) };

    const std::pair<real_type, real_type>
      result
      {
        boost::math::tools::toms748_solve(local::my_riemann_function, bt_val.first, bt_val.second, tol, max_iter)
      };

    {
      std::stringstream strm;

      strm << std::setprecision(std::numeric_limits<real_type>::digits10) << result.first;

      std::cout << strm.str() << std::endl;
    }

    static const real_type my_half { real_type { static_cast<int>(INT8_C(1)) } / static_cast<int>(INT8_C(2)) };

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

  // rz0: 14.134725141734693790457251983562470270784257115699243175685567460149963429809256764949010393171561012
  // rz1: 21.022039638771554992628479593896902777334340524902781754629520403587598586068890799713658514180151419
  // rz2: 25.01085758014568876321379099256282181865954967255799667249654200674509209844164427784023822455806244
  // rz3: 30.424876125859513210311897530584091320181560023715440180962146036993329389333277920290584293902089105
  // rz4: 32.935061587739189690662368964074903488812715603517039009280003440784815608630551005938848496135348718
  // rz5: 37.586178158825671257217763480705332821405597350830793218333001113622149089618537264730329104945823796
  // rz6: 40.918719012147495187398126914633254395726165962777279536161303667253280528720071282996003719889546875

  using local::complex_type;
  using local::find_riemann_root;

  const complex_type rz0 { find_riemann_root(bt_val0) };
  const complex_type rz1 { find_riemann_root(bt_val1) };
  const complex_type rz2 { find_riemann_root(bt_val2) };
  const complex_type rz3 { find_riemann_root(bt_val3) };
  const complex_type rz4 { find_riemann_root(bt_val4) };
  const complex_type rz5 { find_riemann_root(bt_val5) };
  const complex_type rz6 { find_riemann_root(bt_val6) };

  const real_type my_real_tol { std::numeric_limits<real_type>::epsilon() * 256 };

  using std::abs;

  const auto result_rz0_is_ok = (abs(riemann_zeta(rz0)) < my_real_tol);
  const auto result_rz1_is_ok = (abs(riemann_zeta(rz1)) < my_real_tol);
  const auto result_rz2_is_ok = (abs(riemann_zeta(rz2)) < my_real_tol);
  const auto result_rz3_is_ok = (abs(riemann_zeta(rz3)) < my_real_tol);
  const auto result_rz4_is_ok = (abs(riemann_zeta(rz4)) < my_real_tol);
  const auto result_rz5_is_ok = (abs(riemann_zeta(rz5)) < my_real_tol);
  const auto result_rz6_is_ok = (abs(riemann_zeta(rz6)) < my_real_tol);

  const auto
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
