///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2024.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <example/zeta_detail.h>
#include <extended_complex.h>
#include <util.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

auto example023_riemann_zeta_z() -> bool
{
  constexpr unsigned multiprecision_digits10 { static_cast<unsigned>(UINT16_C(101)) };

  using multiprecision_float_type =
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<multiprecision_digits10>,
                                  boost::multiprecision::et_off>;

  using complex_type = extended_complex::complex<multiprecision_float_type>;
  using real_type    = typename complex_type::value_type;

  // N[Zeta[(11/10) + ((23 I) /10), 101]
  //   0.632109498389343535342169571825433072547166503805556530807503296226452975136793086061836360968673834125891496185123905777760
  // - 0.265505793636743413619960696457985203582955058856950038898137949405729351965402359763549645860763401989286376534257444945731 I

  complex_type
    cpx
    {
      real_type(static_cast<int>(INT8_C(11))) / static_cast<int>(INT8_C(10)),
      real_type(static_cast<int>(INT8_C(23))) / static_cast<int>(INT8_C(10))
    };

  const auto rz = riemann_zeta(cpx);

  {
    std::stringstream strm { };

    strm << std::setprecision(std::numeric_limits<real_type>::digits10) << rz;

    std::cout << strm.str() << std::endl;
  }

  bool result_is_ok { };

  {
    const real_type tol { std::numeric_limits<real_type>::epsilon() * static_cast<unsigned>(UINT8_C(64)) };

    const complex_type
      rz_ctrl
      {
        real_type("+0.632109498389343535342169571825433072547166503805556530807503296226452975136793086061836360968673834125891496185123905777760"),
        real_type("-0.265505793636743413619960696457985203582955058856950038898137949405729351965402359763549645860763401989286376534257444945731")
      };

    const auto result_zeta_real_is_ok = util::is_close_fraction(rz.real(), rz_ctrl.real(), tol);
    const auto result_zeta_imag_is_ok = util::is_close_fraction(rz.imag(), rz_ctrl.imag(), tol);

    result_is_ok = (result_zeta_real_is_ok && result_zeta_imag_is_ok);

    std::stringstream strm { };

    strm << "result_is_ok: " << std::boolalpha << result_is_ok;

    std::cout << strm.str() << std::endl;
  }

  return result_is_ok;
}

