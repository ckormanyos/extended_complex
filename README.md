extended_complex
==================

<p align="center">
    <a href="https://github.com/ckormanyos/extended_complex/actions">
        <img src="https://github.com/ckormanyos/extended_complex/actions/workflows/extended_complex.yml/badge.svg" alt="Build Status"></a>
    <a href="https://github.com/ckormanyos/extended_complex/blob/main/LICENSE_1_0.txt">
        <img src="https://img.shields.io/badge/license-BSL%201.0-blue.svg" alt="Boost Software License 1.0"></a>
    <a href="https://godbolt.org/z/E4jz4s43r" alt="godbolt">
        <img src="https://img.shields.io/badge/try%20it%20on-godbolt-green" /></a>
</p>

`ckormanyos/extended_complex` creates an extended-complex number-adaption class.
The project is written in header-only C++14, and compatible through C++14, 17, 20, 23 and beyond..

The `extended_complex::complex` template class can be used
with both built-in floating-point types as well as user-defined numeric types.

## Example

The following straightforward example takes a user-defined,
multiple-precision floating-point type from `Boost.Multiprecision`.
It computes a complex-valued square root with
${\sim}100$ decimal digits of precision.

The square root value computed is

$$
\sqrt { \frac{12}{10} + \frac{34}{10}i }
$$

the approximate complex-value of which is

$$
1.550088912847258141616{\ldots}~{+}~1.096711282759503047577{\ldots}i{\text{.}}
$$

The example code is listed in its entirety below. It is also available _live_
at [Godbolt](https://godbolt.org/z/E4jz4s43r).

```cpp
#include <extended_complex.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <iomanip>
#include <iostream>

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

## In-Depth Example

### Complex-Valued Riemann-Zeta Function

<p align="center">
    <a href="https://godbolt.org/z/scqq9jY1b" alt="godbolt">
        <img src="https://img.shields.io/badge/try%20it%20on-godbolt-green" /></a>
</p>

An in-depth, non-trivial [example](https://github.com/ckormanyos/extended_complex/blob/main/example/example023_riemann_zeta_z.cpp)
provides a header-only implementation of the complex-valued
[Riemann-zeta function](https://github.com/ckormanyos/extended_complex/blob/main/example/zeta_detail.h).
The program handles arguments in a relatively large, yet limited
unit disc of radius ${\sim}~{10}^{6}$ in ${\mathbb{C}}$.
This example uses the algorithm described and found in
the (legacy) `e_float` code and [paper](https://doi.acm.org/10.1145/1916461.1916469).
See also [1] in the references below.

The zeta-function calculation for a single complex-valued point having $101$
decimal digits of precision can also be seen
[here](https://godbolt.org/z/scqq9jY1b).

In particular, the value of

$${\zeta}{\Bigl(}\frac{11}{10} + \frac{23}{10}i{\Bigr)}~{\approx}~0.632109498389343535342{\ldots}~{-}~0.265505793636743413620{\ldots} i$$

is calculated.

### The First Three High-Precision Zeros on the Critical Strip

In [example023a_riemann_zeta_zeros.cpp](https://github.com/ckormanyos/extended_complex/blob/main/example/example023a_riemann_zeta_zeros.cpp),
the first three non-trivial zeros of the complex-valued Riemann-zeta function on the critical strip
are calculated to ${\sim}~101$ decimal digits of precision.

The results found are:

$$
r_{0}~{\sim}~14.13472514173469379045725198356247{\ldots}
$$

$$
r_{1}~{\sim}~21.02203963877155499262847959389690{\ldots}
$$

$$
r_{2}~{\sim}~25.01085758014568876321379099256282{\ldots}
$$

### Not Number-Theory-Ready

The range and domain of the calculations
in this particular example are intended
for high-precision investigations within the above-mentioned unit-disc
of radius ${\sim}~{10}^{6}$ in ${\mathbb{C}}$.

These are _not_ intended for finding record-breaking,
(relatively) low-precision counts of zero-crossings
in the critical strip at ${\mathbb{Re}}(z)=\frac{1}{2}$
which are valuable for providing empirical evidence
for prime number investigations in number-theory.
Other algorithms are needed for this type
of number-theoretical research. See also [2]
for a summary of these and a recent record-breaking
calculation.

## Testing and Continuous Integration

A small test program exercises a variety of non-trivial
algebraic and elementary-function values. The test program verifies
the extended-complex class for both built-in floating point types
`float`, `double` and `long double` as well as a $100$-decimal digit type
from [Boost.Multiprecision](https://www.boost.org/doc/libs/1_84_0/libs/multiprecision/doc/html/index.html).

The above-mentioned in-depth Riemann-zeta example is also executed
and verified in CI.

Continuous integration runs on Ubuntu with both GCC/clang using
the develop branch of modular-boost.

## References

[1] C.M. Kormanyos,
_Algorithm_ _910_: _A_ _Portable_ _C++_ _Multiple_-_Precision_ _System_ _for_ _Special_-_Function_ _Calculations_,
ACM Transactions on Mathematical Software, Vol. 37, Issue 4, pp 1-27.
See also the following [link](https://doi.org/10.1145/1916461.1916469).

[2] D. Platt and T. Trudgian,
_The_ _Riemann_ _hypothesis_ _is_ _true_ _up_ _to_ ${\mbox{\textit{3}}}{\cdot}{\mbox{\textit{10}}}^{\mbox{\textit{\small{12}}}}$,
[arXiv:2004.09765](https://arxiv.org/pdf/2004.09765.pdf) [math.NT].
