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
#if defined(EXTENDED_COMPLEX_RIEMANN_USE_STD_COMPLEX)
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
#include <complex>
#endif
#include <iomanip>
#include <iostream>
#include <sstream>

auto example023_riemann_zeta_z() -> bool
{
  constexpr unsigned multiprecision_digits10 { static_cast<unsigned>(UINT16_C(101)) };

  using multiprecision_float_type =
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<multiprecision_digits10>,
                                  boost::multiprecision::et_off>;

  #if !defined(EXTENDED_COMPLEX_RIEMANN_USE_STD_COMPLEX)
  using complex_type = extended_complex::complex<multiprecision_float_type>;
  #else
  using complex_type = std::complex<multiprecision_float_type>;
  #endif
  using real_type    = typename complex_type::value_type;

  // N[Zeta[(11/10) + ((23 I) /10), 101]
  //   0.632109498389343535342169571825433072547166503805556530807503296226452975136793086061836360968673834125891496185123905777760
  // - 0.265505793636743413619960696457985203582955058856950038898137949405729351965402359763549645860763401989286376534257444945731 I

  complex_type c(real_type(11) / 10, real_type(23) / 10);

  const auto rz = riemann_zeta(c);

  {
    std::stringstream strm;

    strm << std::setprecision(std::numeric_limits<real_type>::digits10) << rz;

    std::cout << strm.str() << std::endl;
  }

  bool result_is_ok { };

  {
    const real_type tol { std::numeric_limits<real_type>::epsilon() * 64U };

    const complex_type
      rz_ctrl
      {
        real_type("+0.632109498389343535342169571825433072547166503805556530807503296226452975136793086061836360968673834125891496185123905777760"),
        real_type("-0.265505793636743413619960696457985203582955058856950038898137949405729351965402359763549645860763401989286376534257444945731")
      };

    const auto result_zeta_real_is_ok = util::is_close_fraction(rz.real(), rz_ctrl.real(), tol);
    const auto result_zeta_imag_is_ok = util::is_close_fraction(rz.imag(), rz_ctrl.imag(), tol);

    result_is_ok = (result_zeta_real_is_ok && result_zeta_imag_is_ok);

    std::stringstream strm;

    strm << "result_is_ok: " << std::boolalpha << result_is_ok;

    std::cout << strm.str() << std::endl;
  }

  return result_is_ok;
}

#if defined(EXTENDED_COMPLEX_RIEMANN_USE_STD_COMPLEX)
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
#endif
