extended_complex
==================

<p align="center">
    <a href="https://github.com/ckormanyos/extended_complex/actions">
        <img src="https://github.com/ckormanyos/extended_complex/actions/workflows/extended_complex.yml/badge.svg" alt="Build Status"></a>
    <a href="https://github.com/ckormanyos/extended_complex/blob/main/LICENSE_1_0.txt">
        <img src="https://img.shields.io/badge/license-BSL%201.0-blue.svg" alt="Boost Software License 1.0"></a>
    <a href="https://godbolt.org/z/hrvhTsGas" alt="godbolt">
        <img src="https://img.shields.io/badge/try%20it%20on-godbolt-green" /></a>
</p>

`extended_complex` creates an extended-complex-number adaption-class in C++14.

It can be used with both built-in floating-point types as well as
user-defined numberic types.

## Example

The following example uses a user-defined multiple-precision floating-point
type from `Boost.Multiprecision` to compute the value of a square root
with approximately $100$ decinal digits of precision

The approximate square root value computed is:

$$
\sqrt {\left(} 1.2+3.4i {\right)} {\sim}1.550088912847258141616125654603881566976{\ldots}+1.096711282759503047577277387056220643003{\ldots}i{\text{.}}
$$

The example code is listed in its entirety below. It is also available _live_
at [Godbolt](https://godbolt.org/z/hrvhTsGas).

```cpp
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include <extended_complex.h>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace local
{
  template<typename NumericType>
  auto is_close_fraction(const NumericType& a,
                         const NumericType& b,
                         const NumericType& tol = std::numeric_limits<NumericType>::epsilon() * 64) noexcept -> bool
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
} // namespace local

auto main() -> int
{
  using float_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>, boost::multiprecision::et_off>;

  using complex_type = extended_complex::complex<float_type>;

  const complex_type val_z1(float_type(12U) / 10U, float_type(34U) / 10U);

  const complex_type sqrt_result = sqrt(val_z1);

  const char* p_str_ctrl_real = "+1.5500889128472581416161256546038815669761567486848749301860666965618993040312647033986371788677357208";
  const char* p_str_ctrl_imag = "+1.096711282759503047577277387056220643003106823143745046422869808875853261131777962620301480493467395";

  const auto result_real_is_ok = local::is_close_fraction(sqrt_result.real(), float_type(p_str_ctrl_real));
  const auto result_imag_is_ok = local::is_close_fraction(sqrt_result.imag(), float_type(p_str_ctrl_imag));

  const auto result_is_ok = (result_real_is_ok && result_imag_is_ok);

  // Print the hexadecimal representation string output.
  const auto flg = std::cout.flags();

  // Visualize if the result is OK.
  std::cout << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10))
            << sqrt_result
            << std::endl;

  // Print the exact result-OK indication.
  std::cout << "result_is_ok: " << std::boolalpha << result_is_ok << std::endl;

  std::cout.flags(flg);
}
```

## Testing and Continuous Integration

A small test program exercises a variety of non-trivial
algebraic and elementary-function values. The test program tests
both built-in floating point types `float`, `double`
and `long double` as well as a $100$-decimal digit type from
`Boost.Multiprecision`.

Continuous integration runs on Ubuntu with both GCC/clang using
the develop branch of modular-boost.
