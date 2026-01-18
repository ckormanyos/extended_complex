///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2024 - 2026.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <example/zeta_detail.h>
#include <extended_complex.h>
#include <util.h>

namespace local
{
  auto riemann_zeta_z_check_double() -> bool;

  auto riemann_zeta_z_check_double() -> bool
  {
    using complex_type = extended_complex::complex<double>;
    using real_type    = typename complex_type::value_type;

    // N[Zeta[(11/10) + ((23 I) /10), 101]
    //   0.632109498389343535342169571825433072547166503805556530807503296226452975136793086061836360968673834125891496185123905777760
    // - 0.265505793636743413619960696457985203582955058856950038898137949405729351965402359763549645860763401989286376534257444945731 I

    const complex_type
      cpx
      {
        real_type(static_cast<int>(INT8_C(11))) / static_cast<int>(INT8_C(10)),
        real_type(static_cast<int>(INT8_C(23))) / static_cast<int>(INT8_C(10))
      };

    const auto rz = riemann_zeta(cpx);

    bool result_is_ok { };

    {
      const real_type tol { std::numeric_limits<real_type>::epsilon() * static_cast<unsigned>(UINT8_C(64)) };

      const complex_type
        rz_ctrl
        {
          real_type { +0.6321094983893435353421695718254330725471665038055565 },
          real_type { -0.2655057936367434136199606964579852035829550588569500 }
        };

      const auto result_zeta_real_is_ok = util::is_close_fraction(rz.real(), rz_ctrl.real(), tol);
      const auto result_zeta_imag_is_ok = util::is_close_fraction(rz.imag(), rz_ctrl.imag(), tol);

      result_is_ok = (result_zeta_real_is_ok && result_zeta_imag_is_ok);
    }

    return result_is_ok;
  }
} // namespace local

auto example023b_riemann_zeta_z_check_double() -> bool
{
  const bool result_check_double_is_ok { local::riemann_zeta_z_check_double() };

  return result_check_double_is_ok;
}
