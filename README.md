extended_complex
==================

<p align="center">
    <a href="https://github.com/ckormanyos/extended_complex/actions">
        <img src="https://github.com/ckormanyos/extended_complex/actions/workflows/extended_complex.yml/badge.svg" alt="Build Status"></a>
    <a href="https://github.com/ckormanyos/extended_complex/blob/main/LICENSE_1_0.txt">
        <img src="https://img.shields.io/badge/license-BSL%201.0-blue.svg" alt="Boost Software License 1.0"></a>
    <a href="https://godbolt.org/z/Paadqoqbv" alt="godbolt">
        <img src="https://img.shields.io/badge/try%20it%20on-godbolt-green" /></a>
</p>

`ckormanyos/extended_complex` creates an extended complex-number adaption class.
The project is written in header-only C++14, and compatible through C++14, 17, 20, 23 and beyond.

The `extended_complex::complex` template class can be used
with both built-in floating-point types as well as user-defined numeric types.

## Example

The following straightforward example takes a user-defined,
multiple-precision floating-point type from
[Boost.Multiprecision](https://www.boost.org/doc/libs/1_84_0/libs/multiprecision/doc/html/index.html).
It computes a complex-valued square root with
${\sim}~100$ decimal digits of precision.

The square root value computed is

$$
    \sqrt { \frac{12}{10} + \frac{34}{10}i }
$$

$$
{\approx}~1.550088912847258141616{\ldots}~{+}~1.096711282759503047577{\ldots}i{\text{.}}
$$

The example code is listed in its entirety below. It is also available _live_
at [Godbolt](https://godbolt.org/z/Paadqoqbv).

```cpp
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <extended_complex.h>

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

    return (fabs(1 - (a / b)) < tol);
  }
} // namespace local

auto main() -> int
{
  using complex_type = extended_complex::complex<boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>, boost::multiprecision::et_off>>;

  using real_type = typename complex_type::value_type;

  const complex_type val_z1(real_type(12U) / 10U, real_type(34U) / 10U);

  const complex_type sqrt_result { sqrt(val_z1) };

  const real_type ctrl_real { "+1.5500889128472581416161256546038815669761567486848749301860666965618993040312647033986371788677357208" };
  const real_type ctrl_imag { "+1.096711282759503047577277387056220643003106823143745046422869808875853261131777962620301480493467395" };

  const auto result_is_ok = (   local::is_close_fraction(sqrt_result.real(), ctrl_real)
                             && local::is_close_fraction(sqrt_result.imag(), ctrl_imag));

  // Print the hexadecimal representation string output.
  const auto flg = std::cout.flags();

  // Visualize if the result is OK.
  std::cout << std::setprecision(static_cast<std::streamsize>(std::numeric_limits<real_type>::digits10))
            << sqrt_result
            << std::endl;

  // Print the result-OK indication.
  std::cout << "result_is_ok: " << std::boolalpha << result_is_ok << std::endl;

  std::cout.flags(flg);
}
```

## In-Depth Example

### Complex-Valued Riemann-Zeta Function

<p align="center">
    <a href="https://godbolt.org/z/xqcehj4fj" alt="godbolt">
        <img src="https://img.shields.io/badge/try%20it%20on-godbolt-green" /></a>
</p>

An in-depth, non-trivial [example](https://github.com/ckormanyos/extended_complex/blob/main/example/example023_riemann_zeta_z.cpp)
provides a header-only implementation of the complex-valued
[Riemann-zeta function](https://github.com/ckormanyos/extended_complex/blob/main/example/zeta_detail.h).
The program handles arguments in a relatively large, yet limited,
unit disc of radius ${\sim}~{10}^{6}$ in ${\mathbb{C}}$.
This example uses the algorithm described and found in
the [`e_float`](https://doi.acm.org/10.1145/1916461.1916469)
code and paper. See also [1] in the references below.

The zeta-function calculation for a single complex-valued point having $101$
decimal digits of precision can also be seen
[here](https://godbolt.org/z/xqcehj4fj).

In particular, the value of

$$
{\zeta}{\Bigl(}\frac{11}{10} + \frac{23}{10}i{\Bigr)}
$$

$$
{\approx}{~}{~}0.632109498389343535342{\ldots}~{-}~0.265505793636743413620{\ldots} i$$
$$

is calculated.

### High-Precision Zeros on the Critical Line

The so-called _critical_ _strip_ refers to the region in the complex plane
for which the argument $z$ of the complex-valued Riemann-zeta
function ${\zeta}(z)$ is

$$
z={\sigma}~+~it{\mbox{,}}
$$

where $0<{\sigma}<1$.

It is believed that there are infinitely many non-trivial roots (zeros)
of the complex-valued Riemann-zeta function. The critical strip
is thought to contain all the non-trivial zeros.

This characteristic of the Riemann-zeta function forges deep connections
to both prime numbers as well as number theory.

The [Riemann Hypothesis](https://en.wikipedia.org/wiki/Riemann_hypothesis),
for instance, states that all non-trivial zeros
of the Riemann zeta function are further localized and lie on the _critical_ _line_
at ${\sigma}~=~1/2$. This bold conjecture remains unproven despite significant efforts
by mathematicians to prove it.

The graph below shows the absolute value of the complex-valued
Riemann-zeta function on a small segment of the critical line.
The first $7$ non-trivial zeros are visible.

![](./images/zeta_critical_strip.jpg)

This image showing $\vert\zeta(z)\vert$ on a small segment of
the critical line has been obtained from
[WolframAlpha(R)](https://www.wolframalpha.com/input?i=Plot%5BAbs%5BZeta%5B%281%2F2%29+%2B+%28I+t%29%5D%5D%2C+%7Bt%2C+1%2C+42%7D%5D)
using the following command.

```wl
Plot[Abs[Zeta[(1/2) + (I t)]], {t, 1, 42}]
```

In [example023a_riemann_zeta_zeros.cpp](https://github.com/ckormanyos/extended_complex/blob/main/example/example023a_riemann_zeta_zeros.cpp),
the first $7$ non-trivial zeros of the complex-valued Riemann-zeta function on the critical line
are calculated to ${\sim}~501$ decimal digits of precision.
Root finding uses [_Algorithm_ _748_](https://doi.org/10.1145/210089.210111).
See also [2] in the references below.

The results found are:

$$
t_{0}~{\approx}~14.134725141734693790457251983562470270784257115699{\ldots}
$$

$$
t_{1}~{\approx}~21.022039638771554992628479593896902777334340524902{\ldots}
$$

$$
t_{2}~{\approx}~25.010857580145688763213790992562821818659549672557{\ldots}
$$

$$
t_{3}~{\approx}~30.424876125859513210311897530584091320181560023715{\ldots}
$$

$$
t_{4}~{\approx}~32.935061587739189690662368964074903488812715603517{\ldots}
$$

$$
t_{5}~{\approx}~37.586178158825671257217763480705332821405597350830{\ldots}
$$

$$
t_{6}~{\approx}~40.918719012147495187398126914633254395726165962777{\ldots}
$$

### _Not_ Number-Theory-Ready

The range and domain of the Riemann-zeta calculations
in the particular examples described above are designed
for high-precision investigations within the above-mentioned unit-disc
of radius ${\sim}~{10}^{6}$ in ${\mathbb{C}}$.

These are _not_ intended for finding record-breaking,
relatively low-precision counts of zero-crossings
on the critical line at ${\sigma}=1/2$.
These are valuable for providing empirical evidence
for prime number investigations in number-theory.
Other algorithms are needed for this type
of number-theoretical research. See also [3]
for a summary of these methods and a recent
record-breaking calculation.

## Testing

A small test program exercises a variety of non-trivial
algebraic and elementary-function values. The test program verifies
the extended-complex class for both built-in floating point types
`float`, `double` and `long double` as well as a $100$-decimal digit type
from [Boost.Multiprecision](https://www.boost.org/doc/libs/1_84_0/libs/multiprecision/doc/html/index.html).

The above-mentioned in-depth Riemann-zeta examples are also executed
and verified in the tests.

## Continuous Integration

Continuous integration runs on Ubuntu and MacOS with both GCC/clang
and also runs on Windows with MSVC. GCC's run-time
[sanitizers](https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html)
are also used in CI in order to help assure dynamic quality.
CI uses the develop branch of modular-boost, when needed, for multiprecision types.

## Additonal details

`ckormanyos/extended_complex` has been tested with numerous compilers,
including target systems ranging from eight to sixty-four bits.
The library is specifically designed for dualistic efficiency and portability.

### Configuration macros (compile-time)

Various configuration features can optionally be
enabled or disabled at compile time with the compiler switches:

```cpp
#define EXTENDED_COMPLEX_DISABLE_IOSTREAM
```

When working with even the most tiny microcontroller systems,
I/O streaming can optionally be disabled with the compiler switch:

```cpp
#define EXTENDED_COMPLEX_DISABLE_IOSTREAM
```

The default setting is `EXTENDED_COMPLEX_DISABLE_IOSTREAM` not set
and I/O streaming operations are enabled.

```cpp
#define EXTENDED_COMPLEX_CONSTEXPR
```

The macro `EXTENDED_COMPLEX_CONSTEXPR` is default-defined to be equal
to the word `constexpr`. This macro was previously used (for old compilers
no longer supported) to either use or un-use the word `constexpr`.
This was back when `constexpr` was new. At this time, simply leave this
macro unchanged and equal to the word `constexpr`.

```cpp
#define EXTENDED_COMPLEX_RIEMANN_USE_STD_COMPLEX
```

Define this advanced design macro on the command line in order to use
`std::complex` instead of the extended-complex implementation of complex
when performing stress-tests with Riemann-Zeta calculations.
This macro can be used to help ensure that this library is compatible to
the standard-library's `std::complex`.

## References

[1] C.M. Kormanyos,
_Algorithm_ _910_: _A_ _Portable_ _C++_ _Multiple_-_Precision_ _System_ _for_ _Special_-_Function_ _Calculations_,
ACM Transactions on Mathematical Software, Vol. 37, Issue 4, pp 1-27 (01 February 2011),
[https://doi.org/10.1145/1916461.1916469](https://doi.org/10.1145/1916461.1916469).

[2] G.E. Alefeld, F.A. Potra, Yixun Shi,
_Algorithm_ _748_: _enclosing_ _zeros_ _of_ _continuous_ _functions_,
ACM Transactions on Mathematical Software, Vol. 21, Issue 3, pp 327-344 (01 September 1995),
[https://doi.org/10.1145/210089.210111](https://doi.org/10.1145/210089.210111).

[3] D. Platt and T. Trudgian,
_The_ _Riemann_ _hypothesis_ _is_ _true_ _up_ _to_ ${\mbox{\textit{3}}}{\cdot}{\mbox{\textit{10}}}^{\mbox{\textit{\small{12}}}}$,
[arXiv:2004.09765](https://arxiv.org/pdf/2004.09765.pdf).
