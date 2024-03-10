///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2024.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <vector>

#define BOOST_MATH_STANDALONE
#define BOOST_MULTIPRECISION_STANDALONE

#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <extended_complex.h>

namespace local { namespace detail {

namespace Util {

template<typename T1,
         typename T2 = T1>
struct point
{
  T1 my_x { };
  T2 my_y { };

  explicit point(const T1& x = T1(),
                 const T2& y = T2()) : my_x(x), my_y(y) { }
};

} // namespace Util

namespace ef {

void prime_factors(const std::uint32_t n, std::deque<Util::point<std::uint32_t> >& pf);
void prime(const std::uint32_t n, std::deque<std::uint32_t>& primes);

template<typename T>
int order(const T& val)
{
  using std::log;
  using std::lround;

  const auto d10_scale =
    static_cast<std::uint32_t>
    (
      lround
      (
        static_cast<float>
        (
            static_cast<float>
            (
              1000.0F * log(static_cast<float>(std::numeric_limits<T>::radix))
            )
          / log(10.0F)
        )
      )
    );

  const auto ib = static_cast<std::int32_t>(ilogb(val));

  const auto my_order =
    static_cast<std::uint32_t>
    (
        static_cast<std::uint64_t>(static_cast<std::uint64_t>(ib) * d10_scale)
      / static_cast<std::uint32_t>(UINT16_C(1000))
    );

  return my_order;
}

template<typename T>
int tol(void)
{
  return std::numeric_limits<T>::max_digits10;
}

template<typename T>
constexpr T one()
{
  return T { 1U };
}

template<typename T>
constexpr double to_double(const T& val)
{
  return static_cast<double>(val);
}

template<typename T>
constexpr double to_int64(const T& val)
{
  return static_cast<double>(val);
}

template<typename T>
constexpr T int64_min()
{
  return T { (std::numeric_limits<std::int64_t>::min)() };
}

template<typename T>
constexpr T int64_max()
{
  return T { (std::numeric_limits<std::int64_t>::max)() };
}

} // namespace ef

namespace Util {

namespace detail {

template<typename T>
struct logn_helper
{
public:
  static T my_logn(const std::uint32_t n)
  {
    const auto it_ln =
      std::find_if
      (
        ln_data.cbegin(),
        ln_data.cend(),
        [&n](const auto& elem)
        {
          return (elem.first == n);
        }
      );

    if(it_ln == ln_data.cend())
    {
      using std::log;

      const T ln_value = log( T { n } );

      ln_data[n] = ln_value;

      return ln_value;
    }
    else
    {
      return it_ln->second;
    }
  }

private:
  static std::map<std::uint32_t, T> ln_data;
};

template<typename T>
std::map<std::uint32_t, T> logn_helper<T>::ln_data;

} // namespace detail

template<typename T>
T logn(const std::uint32_t n)
{
  return detail::logn_helper<T>::my_logn(n);
}

template<typename ComplexType>
ComplexType j_pow_x(const std::uint32_t j, const ComplexType& x, std::map<std::uint32_t, ComplexType>& n_pow_x_prime_factor_map)
{
  using local_complex_type = ComplexType;
  using local_real_type    = typename local_complex_type::value_type;

  std::deque<Util::point<std::uint32_t>> pf;

  ef::prime_factors(j, pf);

  local_complex_type jpx = ef::one<local_real_type>();

  for(std::size_t i = static_cast<std::size_t>(0u); i < pf.size(); i++)
  {
    local_complex_type pf_pow_x;

    const std::uint32_t n = pf[i].my_x;
    const std::uint32_t p = pf[i].my_y;

    const typename std::map<std::uint32_t, local_complex_type>::const_iterator it = n_pow_x_prime_factor_map.find(n);

    if(it == n_pow_x_prime_factor_map.end())
    {
      // Compute n^x using exp[x * log(n)] and use the map data in the Zeta::logn(...).
      // Obtain the necessary integer logarithms from a table.

      const auto x_is_int = ((x.imag() == 0) && (static_cast<int>(x.real()) == x.real()));

      if(x_is_int)
      {
        const local_real_type rx = x.real();

        // Compute pure integer power for pure integer arguments.
        if((rx < ef::int64_max<local_real_type>()) && (rx > ef::int64_min<local_real_type>()))
        {
          pf_pow_x = pow(local_complex_type(n), static_cast<int>(ef::to_int64(rx)));
        }
        else
        {
          pf_pow_x = exp(x * logn<local_real_type>(n));
        }
      }
      else
      {
        pf_pow_x = exp(x * logn<local_real_type>(n));
      }

      n_pow_x_prime_factor_map[n] = pf_pow_x;
    }
    else
    {
      pf_pow_x = it->second;
    }

    // Do the power expansion.
    if     (p == static_cast<std::uint32_t>(1u)) { }
    else if(p == static_cast<std::uint32_t>(2u)) { pf_pow_x *=  pf_pow_x; }
    else if(p == static_cast<std::uint32_t>(3u)) { pf_pow_x *= (pf_pow_x * pf_pow_x); }
    else                                         { pf_pow_x *= pow(pf_pow_x, static_cast<std::int64_t>(p - 1u)); }

    jpx *= pf_pow_x;
  }

  return jpx;
}

} // namespace Util

namespace Primes {

struct Inserter
{
public:
  static const std::size_t start_index = static_cast<std::size_t>(2u);

  explicit Inserter(std::deque<std::uint32_t>& sequence)
    : count(static_cast<std::uint32_t>(start_index)),
      my_it(std::back_inserter(sequence)) { }

  Inserter() = delete;

  void operator()(const bool& bo_is_not_prime) const
  {
    const bool bo_is_prime = !bo_is_not_prime;

    if(bo_is_prime)
    {
      *my_it = count;
    }

    ++count;
  }

private:
  mutable std::uint32_t count { };
  mutable std::back_insert_iterator<std::deque<std::uint32_t>> my_it;
};

void Generator(const std::uint32_t n, std::deque<std::uint32_t>& primes_data)
{
  // Establish the range of the prime number calculation. Use an approximation
  // related to the prime number theorem to obtain the value of the maximum prime
  // number or a minimum of at least 100. Also be sure to limit this range to
  // within the upper limit of std::uint32_t.

  static const std::uint32_t min_hundred = static_cast<std::uint32_t>(100u);
  static const double xmax        = static_cast<double>((std::numeric_limits<std::uint32_t>::max)());

  const std::uint32_t N     = (std::max)(min_hundred, n);
  const double xn           = static_cast<double>(N);
  const double logn         = ::log(xn);
  const double loglogn      = ::log(logn);
  const double top          = xn * (((logn + loglogn) - 1.0) + ((static_cast<double>(1.8) * loglogn) / logn));
  const double xlim         = (std::min)(top, xmax);
  const std::uint32_t nlim  = static_cast<std::uint32_t>(static_cast<std::uint64_t>(xlim));
  const std::uint32_t limit = (std::max)(n, nlim);

  // Use a sieve algorithm to generate a boolean table representation of the primes.

  std::vector<bool> sieve(static_cast<std::size_t>(limit), false);

  std::uint32_t i = static_cast<std::uint32_t>(Primes::Inserter::start_index);
  std::uint32_t i2;

  while((i2 = static_cast<std::uint32_t>(i * i)) < limit)
  {
    if(!sieve[i])
    {
      for(std::uint32_t j = i2; j < limit; j = static_cast<std::uint32_t>(j + i))
      {
        sieve[j] = true;
      }
    }

    ++i;
  }

  // Extract the prime numbers into the data table by inserting them from the sieve.
  primes_data.clear();

  std::for_each(sieve.begin() + Primes::Inserter::start_index,
                sieve.end(),
                Primes::Inserter(primes_data));

  primes_data.resize(static_cast<std::size_t>(n), static_cast<std::uint32_t>(0u));
}

std::deque<std::uint32_t>& Data(void)
{
  // Create a static data table of primes and return a reference to it.
  static std::deque<std::uint32_t> primes;

  if(primes.empty())
  {
    // Select a maximum count of prime numbers to be stored in the data table.
    // This number is selected such that the value of the highest prime will slightly
    // exceed 0x10000 (decimal 65,536). This number is significant because it is
    // the maximum value which needs to be tested while computing the prime factors
    // of unsigned 32-bit integers, as done in the subroutine Factors(...).
    Primes::Generator(static_cast<std::uint32_t>(6550u), primes);
  }

  return primes;
}

bool IsPrimeFactor(std::uint32_t& np, const std::uint32_t p)
{
  const std::uint32_t q = static_cast<std::uint32_t>(np / p);
  const std::uint32_t r = static_cast<std::uint32_t>(np - static_cast<std::uint32_t>(q * p));

  const bool is_prime_factor = (r == static_cast<std::uint32_t>(0u));

  if(is_prime_factor)
  {
    np = q;
  }
    
  return is_prime_factor;
}

void Factors(const std::uint32_t n, std::deque<Util::point<std::uint32_t> >& pf)
{
  // Compute the prime factors of the unsigned integer n. Use the divide algorithm of
  // "The Art of Computer Programming Volume 2 Semi-numerical Algorithms Third Edition",
  // Donald Knuth (Algorithm A, Chapter 4.5.4, page 380 and pages 378-417).
  static const std::size_t sz = Data().size();

  pf.clear();

  const std::uint32_t sqrt_n = static_cast<std::uint32_t>(static_cast<std::uint64_t>(::sqrt(static_cast<double>(n)) + 0.5));

  std::uint32_t np = n;

  for(std::size_t i = static_cast<std::size_t>(0u); i < sz; i++)
  {
    const std::uint32_t p = Data()[i];

    if(IsPrimeFactor(np, p))
    {
      Util::point<std::uint32_t> ip(p, static_cast<std::uint32_t>(1u));

      while(IsPrimeFactor(np, p))
      {
        ++ip.my_y;
      }

      pf.push_back(ip);
    }

    if(static_cast<std::uint32_t>(np / p) <= p)
    {
      pf.push_back(Util::point<std::uint32_t>(np, static_cast<std::uint32_t>(1u)));

      break;
    }

    if((np == static_cast<std::uint32_t>(1u)) || (p >= sqrt_n))
    {
      break;
    }
  }
}

} // namespace Primes

namespace ef {

void prime(const std::uint32_t n, std::deque<std::uint32_t>& primes)
{
  // For small values of n less than the size of the prime data table, the primes
  // can be copied from the data table. For large values of n, the primes must be
  // generated.
  if(n < static_cast<std::uint32_t>(Primes::Data().size()))
  {
    primes.assign(Primes::Data().begin(), Primes::Data().begin() + static_cast<std::size_t>(n));
  }
  else
  {
    Primes::Generator(n, primes);
  }
}

void prime_factors(const std::uint32_t n, std::deque<Util::point<std::uint32_t>>& pf)
{
  // Factor the input integer into a list of primes. For small inputs less than 10,000
  // use the tabulated prime factors list. Calculate the prime factors for larger inputs
  // above 10,000.
  static std::vector<std::deque<Util::point<std::uint32_t> > > prime_factors_list;

  if(prime_factors_list.empty())
  {
    // Generate a table of the sets of the first 10,000 integer prime factorizations.
    prime_factors_list.resize(static_cast<std::size_t>(10000u));

    prime_factors_list[static_cast<std::size_t>(0u)] = std::deque<Util::point<std::uint32_t> >(static_cast<std::size_t>(1u), Util::point<std::uint32_t>(static_cast<std::uint32_t>(0u), static_cast<std::uint32_t>(1u)));
    prime_factors_list[static_cast<std::size_t>(1u)] = std::deque<Util::point<std::uint32_t> >(static_cast<std::size_t>(1u), Util::point<std::uint32_t>(static_cast<std::uint32_t>(1u), static_cast<std::uint32_t>(1u)));
    prime_factors_list[static_cast<std::size_t>(2u)] = std::deque<Util::point<std::uint32_t> >(static_cast<std::size_t>(1u), Util::point<std::uint32_t>(static_cast<std::uint32_t>(2u), static_cast<std::uint32_t>(1u)));
    prime_factors_list[static_cast<std::size_t>(3u)] = std::deque<Util::point<std::uint32_t> >(static_cast<std::size_t>(1u), Util::point<std::uint32_t>(static_cast<std::uint32_t>(3u), static_cast<std::uint32_t>(1u)));

    static const std::uint32_t n_five = static_cast<std::uint32_t>(5u);

    std::deque<std::uint32_t>::const_iterator it_next_prime = std::find(Primes::Data().cbegin(), Primes::Data().cend(), n_five);

    for(std::size_t i = static_cast<std::size_t>(4u); i < prime_factors_list.size(); i++)
    {
      if((it_next_prime != Primes::Data().cend()) && (static_cast<std::uint32_t>(i) == *it_next_prime))
      {
        ++it_next_prime;

        prime_factors_list[i] = std::deque<Util::point<std::uint32_t> >(static_cast<std::size_t>(1u),
                                                                 Util::point<std::uint32_t>(static_cast<std::uint32_t>(i),
                                                                 static_cast<std::uint32_t>(1u)));
      }
      else
      {
        Primes::Factors(static_cast<std::uint32_t>(i), prime_factors_list[i]);
      }
    }
  }

  if(static_cast<std::size_t>(n) < prime_factors_list.size())
  {
    pf = prime_factors_list[static_cast<std::size_t>(n)];
  }
  else
  {
    Primes::Factors(n, pf);
  }
}

} // namespace ef

template<typename ComplexType>
ComplexType ZetaTemplate(const ComplexType& s)
{
  using local_complex_type = ComplexType;
  using local_real_type    = typename local_complex_type::value_type;

  //if(ef::isint(s))
  //{
  //  // Support pure-integer arguments according to certain conditions.
  //  const std::int32_t n = ef::to_int32(real(s));
  //
  //  if(Zeta_Series::has_simple_form_for_zeta_n(n))
  //  {
  //    return ef::riemann_zeta(n);
  //  }
  //}

  //if(ef::isneg(s))
  //{
  //  // Support negative arguments using the reflection formula.
  //  // Support arguments with a real part which is less than 1/2.
  //  return Zeta_Series::Reflection(s);
  //}

  // The algorithms for calculating the Riemann zeta function below use calculations
  // of the integer j raised to the power s, or in other words j^s. The calculation of
  // j^s is accelerated using tables of the prime factors of j raised to the power s.
  // The calculation is furthermore accelerated by storing the necessary integer
  // logarithms in a static table.

  // Declare a map of prime factors raised to the power of the argument s.
  std::map<std::uint32_t, local_complex_type> n_pow_s_prime_factor_map;

  // Generate a list of the first 300 prime numbers.
  static std::deque<std::uint32_t> prime_data;

  if(prime_data.empty())
  {
    ef::prime(static_cast<std::uint32_t>(300u), prime_data);
  }

  static const std::vector<std::uint32_t> primes(prime_data.begin(), prime_data.end());

  {
    if(abs(imag(s)) > local_real_type(1000100UL))
    {
      // Return NaN if s has a large imaginary part.
      return local_complex_type { std::numeric_limits<local_real_type>::quiet_NaN() };
    }

    // Use the accelerated alternating converging series for Zeta as shown in:
    // http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.html
    // taken from P. Borwein, "An Efficient Algorithm for the Riemann Zeta Function",
    // January 1995.

    // Compute the coefficients dk in a loop and calculate the zeta function sum
    // within the same loop on the fly.

    // Set up the factorials and powers for the calculation of the coefficients dk.
    // Note that j = n at this stage in the calculation. Also note that the value of
    // dn is equal to the value of d0 at the end of the loop.

    // Use N = (digits * 1.45) + {|imag(s)| * 1.1}
    static const double nd = static_cast<double>(std::numeric_limits<local_real_type>::digits10) * static_cast<double>(1.45);
           const double ni = static_cast<double>(static_cast<double>(1.10) * fabs(ef::to_double(imag(s))));

    const std::int32_t N        = static_cast<std::int32_t>(static_cast<std::int64_t>(static_cast<double>(nd + ni)));
          bool         neg_term = (N % static_cast<std::int32_t>(2)) == static_cast<std::int32_t>(0);

    local_real_type n_plus_j_minus_one_fact = boost::math::factorial<local_real_type>(static_cast<unsigned int>((N + N) - 1));
    local_real_type four_pow_j              = pow(local_real_type(4U), static_cast<std::int64_t>(N));
    local_real_type n_minus_j_fact          = ef::one<local_real_type>();
    local_real_type two_j_fact              = n_plus_j_minus_one_fact * static_cast<std::int32_t>(static_cast<std::int32_t>(2) * N);

    local_real_type dn = (n_plus_j_minus_one_fact * four_pow_j) / (n_minus_j_fact * two_j_fact);

    local_complex_type jps = Util::j_pow_x(static_cast<std::uint32_t>(N), s, n_pow_s_prime_factor_map);

    local_complex_type zs = (!neg_term ? dn : -dn) / jps;

    for(std::int32_t j = N - static_cast<std::int32_t>(1); j >= static_cast<std::int32_t>(0); --j)
    {
      const bool j_is_zero = (j == static_cast<std::int32_t>(0));

      const std::int32_t two_jp1_two_j =
        static_cast<std::int32_t>
        (
            static_cast<std::int32_t>((static_cast<std::int32_t>(2) * j) + static_cast<std::int32_t>(1))
          * static_cast<std::int32_t> (static_cast<std::int32_t>(2) * (!j_is_zero ? j : static_cast<std::int32_t>(1)))
        );

      n_plus_j_minus_one_fact /= static_cast<std::int32_t>(N + j);
      four_pow_j              /= static_cast<std::int32_t>(4);
      n_minus_j_fact          *= static_cast<std::int32_t>(N - j);
      two_j_fact              /= two_jp1_two_j;

      dn += ((n_plus_j_minus_one_fact * four_pow_j) / (n_minus_j_fact * two_j_fact));

      if(!j_is_zero)
      {
        // Increment the zeta function sum.
        jps = Util::j_pow_x(j, s, n_pow_s_prime_factor_map);

        neg_term = !neg_term;

        zs += (!neg_term ? dn : -dn) / jps;
      }
    }

    const local_complex_type two_pow_one_minus_s = pow(local_complex_type( local_real_type { 2 } ), ef::one<local_real_type>() - s);

    return zs / (dn * (ef::one<local_real_type>() - two_pow_one_minus_s));
  }
}

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

} // namespace detail

using complex_type = extended_complex::complex<boost::multiprecision::number<boost::multiprecision::cpp_dec_float<101>, boost::multiprecision::et_off>>;
using real_type    = typename complex_type::value_type;

} // namespace local

template<typename ComplexType>
ComplexType riemann_zeta(const ComplexType& s)
{
  return local::detail::ZetaTemplate<ComplexType>(s);
}

auto example() -> bool
{
  using local::complex_type;
  using local::real_type;

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

    const auto result_zeta_real_is_ok = local::detail::is_close_fraction(rz.real(), rz_ctrl.real(), tol);
    const auto result_zeta_imag_is_ok = local::detail::is_close_fraction(rz.imag(), rz_ctrl.imag(), tol);

    result_is_ok = (result_zeta_real_is_ok && result_zeta_imag_is_ok);

    std::stringstream strm;

    strm << "result_is_ok: " << std::boolalpha << result_is_ok;

    std::cout << strm.str() << std::endl;
  }

  return result_is_ok;
}
