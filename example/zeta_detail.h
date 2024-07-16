///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2024.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef ZETA_DETAIL_2024_03_17_H
#define ZETA_DETAIL_2024_03_17_H

#define BOOST_MATH_STANDALONE
#define BOOST_MULTIPRECISION_STANDALONE

//#define EXTENDED_COMPLEX_RIEMANN_USE_STD_COMPLEX

#include <boost/math/special_functions/factorials.hpp>
#include <boost/unordered/unordered_map.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <vector>

namespace zeta_detail { namespace detail {

namespace Util {

template<typename T1, typename T2 = T1>
struct point
{
  explicit point(const T1& x = T1(),
                 const T2& y = T2()) noexcept : my_x(x), my_y(y) { }

  T1 my_x { };
  T2 my_y { };
};

} // namespace Util

namespace ef {

inline auto prime_factors(const std::uint32_t n, std::deque<Util::point<std::uint32_t> >& pf) -> void;
inline auto prime        (const std::uint32_t n, std::deque<std::uint32_t>& primes)-> void;

template<typename T>
constexpr auto tol(void) noexcept -> int { return std::numeric_limits<T>::max_digits10; }

template<typename T>
constexpr auto one() -> T { return T { 1U }; }

template<typename T>
constexpr auto to_double(const T& val) -> double { return static_cast<double>(val); }

template<typename T>
constexpr auto to_int64(const T& val) -> std::int64_t { return static_cast<std::int64_t>(val); }

template<typename T>
constexpr auto int64_min() -> T { return T { (std::numeric_limits<std::int64_t>::min)() }; }

template<typename T>
constexpr auto int64_max() -> T { return T { (std::numeric_limits<std::int64_t>::max)() }; }

} // namespace ef

namespace Util {

namespace detail {

template<typename T>
struct logn_helper
{
public:
  static auto my_logn(const std::uint32_t n) -> T
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
  static boost::unordered_map<std::uint32_t, T> ln_data;
};

template<typename T>
boost::unordered_map<std::uint32_t, T> logn_helper<T>::ln_data { };

} // namespace detail

template<typename T, typename IntegralType>
auto pown(const T& b, const IntegralType p) -> std::enable_if_t<std::is_integral<IntegralType>::value, T>;

template<typename T>
auto logn(const std::uint32_t n) -> T { return detail::logn_helper<T>::my_logn(n); }

template<typename ComplexType>
auto j_pow_x(const std::uint32_t j, const ComplexType& x, boost::unordered_map<std::uint32_t, ComplexType>& n_pow_x_prime_factor_map) -> ComplexType
{
  using local_complex_type = ComplexType;
  using local_real_type    = typename local_complex_type::value_type;

  std::deque<Util::point<std::uint32_t>> pf;

  ef::prime_factors(j, pf);

  local_complex_type jpx { ef::one<local_real_type>() };

  for(std::size_t i = static_cast<std::size_t>(0u); i < pf.size(); i++)
  {
    local_complex_type pf_pow_x;

    const auto n = pf[i].my_x;
    const auto p = pf[i].my_y;

    using const_iterator_type = typename boost::unordered_map<std::uint32_t, local_complex_type>::const_iterator;

    const const_iterator_type itr = n_pow_x_prime_factor_map.find(n);

    if(itr == n_pow_x_prime_factor_map.end())
    {
      // Compute n^x using exp[x * log(n)] and use the map data in the Zeta::logn(...).
      // Obtain the necessary integer logarithms from a table.

      const bool x_is_int { ((x.imag() == 0) && (static_cast<int>(x.real()) == x.real())) };

      if(x_is_int)
      {
        const local_real_type rx { x.real() };

        // Compute pure integer power for pure integer arguments.
        if((rx < ef::int64_max<local_real_type>()) && (rx > ef::int64_min<local_real_type>()))
        {
          pf_pow_x = Util::pown(local_complex_type(n), ef::to_int64(rx));
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
      pf_pow_x = itr->second;
    }

    // Do the power expansion.
    if     (p == static_cast<std::uint32_t>(UINT8_C(1))) { }
    else if(p == static_cast<std::uint32_t>(UINT8_C(2))) { pf_pow_x *=  pf_pow_x; }
    else if(p == static_cast<std::uint32_t>(UINT8_C(3))) { pf_pow_x *= (pf_pow_x * pf_pow_x); }
    else if(p == static_cast<std::uint32_t>(UINT8_C(4))) { pf_pow_x *= pf_pow_x; pf_pow_x *= pf_pow_x; }
    else                                                 { pf_pow_x *= Util::pown(pf_pow_x, static_cast<std::int64_t>(p - static_cast<std::uint32_t>(UINT8_C(1)))); }

    jpx *= pf_pow_x;
  }

  return jpx;
}

template<typename T, typename IntegralType>
auto pown(const T& b, const IntegralType p) -> std::enable_if_t<std::is_integral<IntegralType>::value, T>
{
  // Calculate (b ^ p) for both real as well as complex-valued b,
  // where p is a signed or unsigned integral type.

  using local_number_type = T;

  local_number_type result { };

  const local_number_type local_one { ef::one<local_number_type>() };

  if     (p <  static_cast<std::int64_t>(INT8_C(0))) { result = local_one / pown(b, -p); }
  else if(p == static_cast<std::int64_t>(INT8_C(0))) { result = local_one; }
  else if(p == static_cast<std::int64_t>(INT8_C(1))) { result = b; }
  else if(p == static_cast<std::int64_t>(INT8_C(2))) { result = b; result *= b; }
  else if(p == static_cast<std::int64_t>(INT8_C(3))) { result = b; result *= b; result *= b; }
  else if(p == static_cast<std::int64_t>(INT8_C(4))) { result = b; result *= b; result *= result; }
  else
  {
    result = local_one;

    local_number_type y(b);

    auto p_local = static_cast<std::uint64_t>(p);

    // Use the so-called ladder method for the power calculation.
    for(;;)
    {
      const auto do_power_multiply =
        (static_cast<std::uint_fast8_t>(p_local & static_cast<unsigned>(UINT8_C(1))) != static_cast<std::uint_fast8_t>(UINT8_C(0)));

      if(do_power_multiply)
      {
        result *= y;
      }

      p_local >>= static_cast<unsigned>(UINT8_C(1));

      if(p_local == static_cast<std::uint64_t>(UINT8_C(0)))
      {
        break;
      }

      y *= y;
    }
  }

  return result;
}

} // namespace Util

namespace Primes {

template<typename IntegralType>
struct Inserter
{
public:
  using value_type = IntegralType;

  static constexpr std::size_t start_index = static_cast<std::size_t>(UINT8_C(2));

  explicit Inserter(std::deque<value_type>& sequence)
    : my_it(std::back_inserter(sequence)) { }

  Inserter() = delete;

  void operator()(const bool& bo_is_not_prime) const
  {
    const auto bo_is_prime = (!bo_is_not_prime);

    if(bo_is_prime)
    {
      *my_it = count;
    }

    ++count;
  }

private:
  mutable value_type count { start_index };
  mutable std::back_insert_iterator<std::deque<value_type>> my_it;
};

inline auto Generator(const std::uint32_t n, std::deque<std::uint32_t>& primes_data) -> void
{
  // Establish the range of the prime number calculation. Use an approximation
  // related to the prime number theorem to obtain the value of the maximum prime
  // number or a minimum of at least 100. Also be sure to limit this range to
  // within the upper limit of std::uint32_t.

  constexpr auto min_hundred = static_cast<std::uint32_t>(UINT8_C(100));
  constexpr auto xmax        = static_cast<double>((std::numeric_limits<std::uint32_t>::max)());

  const auto xn = static_cast<double>((std::max)(min_hundred, n));

  using std::log;

  const auto logn    = log(xn);
  const auto loglogn = log(logn);
  const auto top     = static_cast<double>(xn * (((logn + loglogn) - 1.0) + ((1.8 * loglogn) / logn)));
  const auto xlim    = (std::min)(top, xmax);
  const auto nlim    = static_cast<std::uint32_t>(static_cast<std::uint64_t>(xlim));
  const auto limit   = (std::max)(n, nlim);

  // Use a sieve algorithm to generate a boolean table representation of the primes.

  std::vector<bool> sieve(static_cast<std::size_t>(limit), false);

  using primes_inserter_type = Primes::Inserter<std::uint32_t>;

  auto i = static_cast<std::uint32_t>(primes_inserter_type::start_index);

  std::uint32_t i2 { };

  while((i2 = static_cast<std::uint32_t>(i * i)) < limit)
  {
    if(!sieve[i])
    {
      for(auto j = i2; j < limit; j = static_cast<std::uint32_t>(j + i))
      {
        sieve[j] = true;
      }
    }

    ++i;
  }

  // Extract the prime numbers into the data table by inserting them from the sieve.
  primes_data.clear();

  std::for_each(sieve.begin() + primes_inserter_type::start_index,
                sieve.end(),
                primes_inserter_type(primes_data));

  primes_data.resize(static_cast<std::size_t>(n), static_cast<std::uint32_t>(0u));
}

inline auto Data(void) -> std::deque<std::uint32_t>&
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
    Primes::Generator(static_cast<std::uint32_t>(UINT16_C(6550)), primes);
  }

  return primes;
}

constexpr auto IsPrimeFactor(std::uint32_t& np, const std::uint32_t p) -> bool
{
  const std::uint32_t q = static_cast<std::uint32_t>(np / p);
  const std::uint32_t r = static_cast<std::uint32_t>(np - static_cast<std::uint32_t>(q * p));

  const bool is_prime_factor = (r == static_cast<std::uint32_t>(UINT8_C(0)));

  if(is_prime_factor)
  {
    np = q;
  }

  return is_prime_factor;
}

inline auto Factors(const std::uint32_t n, std::deque<Util::point<std::uint32_t> >& pf) -> void
{
  // Compute the prime factors of the unsigned integer n. Use the divide algorithm of
  // "The Art of Computer Programming Volume 2 Semi-numerical Algorithms Third Edition",
  // Donald Knuth (Algorithm A, Chapter 4.5.4, page 380 and pages 378-417).
  static const std::size_t sz = Data().size();

  pf.clear();

  using std::sqrt;

  const auto sqrt_n = static_cast<std::uint32_t>(static_cast<std::uint64_t>(sqrt(static_cast<double>(n)) + 0.5));

  auto np = n;

  for(auto i = static_cast<std::size_t>(UINT8_C(0)); i < sz; ++i)
  {
    const auto p = Data()[i];

    if(IsPrimeFactor(np, p))
    {
      Util::point<std::uint32_t> ip(p, static_cast<std::uint32_t>(UINT8_C(1)));

      while(IsPrimeFactor(np, p))
      {
        ++ip.my_y;
      }

      pf.push_back(ip);
    }

    if(static_cast<std::uint32_t>(np / p) <= p)
    {
      pf.push_back(Util::point<std::uint32_t>(np, static_cast<std::uint32_t>(UINT8_C(1))));

      break;
    }

    if((np == static_cast<std::uint32_t>(UINT8_C(1))) || (p >= sqrt_n))
    {
      break;
    }
  }
}

} // namespace Primes

namespace ef {

inline auto prime(const std::uint32_t n, std::deque<std::uint32_t>& primes) -> void
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

inline auto prime_factors(const std::uint32_t n, std::deque<Util::point<std::uint32_t>>& pf) -> void
{
  using local_point_type       = Util::point<std::uint32_t>;
  using local_point_deque_type = std::deque<local_point_type>;

  // Factor the input integer into a list of primes. For small inputs less than 10,000
  // use the tabulated prime factors list. Calculate the prime factors for larger inputs
  // above 10,000.
  static std::vector<local_point_deque_type> prime_factors_list { };

  if(prime_factors_list.empty())
  {
    // Generate a table of the sets of the first 10,000 integer prime factorizations.
    prime_factors_list.resize(static_cast<std::size_t>(UINT16_C(10000)));

    prime_factors_list[static_cast<std::size_t>(UINT8_C(0))] = local_point_deque_type(static_cast<std::size_t>(UINT8_C(1)), local_point_type(static_cast<std::uint32_t>(UINT8_C(0)), static_cast<std::uint32_t>(UINT8_C(1))));
    prime_factors_list[static_cast<std::size_t>(UINT8_C(1))] = local_point_deque_type(static_cast<std::size_t>(UINT8_C(1)), local_point_type(static_cast<std::uint32_t>(UINT8_C(1)), static_cast<std::uint32_t>(UINT8_C(1))));
    prime_factors_list[static_cast<std::size_t>(UINT8_C(2))] = local_point_deque_type(static_cast<std::size_t>(UINT8_C(1)), local_point_type(static_cast<std::uint32_t>(UINT8_C(2)), static_cast<std::uint32_t>(UINT8_C(1))));
    prime_factors_list[static_cast<std::size_t>(UINT8_C(3))] = local_point_deque_type(static_cast<std::size_t>(UINT8_C(1)), local_point_type(static_cast<std::uint32_t>(UINT8_C(3)), static_cast<std::uint32_t>(UINT8_C(1))));

    constexpr std::uint32_t n_five { static_cast<std::uint32_t>(UINT8_C(5)) };

    auto it_next_prime = std::find(Primes::Data().cbegin(), Primes::Data().cend(), n_five);

    for(auto i = static_cast<std::size_t>(UINT8_C(4)); i < prime_factors_list.size(); ++i)
    {
      if((it_next_prime != Primes::Data().cend()) && (static_cast<std::uint32_t>(i) == *it_next_prime))
      {
        ++it_next_prime;

        const local_point_type pt_i_one { static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(UINT8_C(1)) };

        prime_factors_list[i] = local_point_deque_type(static_cast<std::size_t>(UINT8_C(1)), pt_i_one);
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
auto ZetaTemplate(const ComplexType& s) -> ComplexType
{
  using local_complex_type = ComplexType;
  using local_real_type    = typename local_complex_type::value_type;

  // TODO: Support pure-integer arguments according to certain conditions.
  // TODO: Support negative arguments using the reflection formula.
  // TODO: Support arguments with a real part which is less than 1/2.

  // The algorithms for calculating the Riemann zeta function below use calculations
  // of the integer j raised to the power s, or in other words j^s. The calculation of
  // j^s is accelerated using tables of the prime factors of j raised to the power s.
  // The calculation is furthermore accelerated by storing the necessary integer
  // logarithms in a static table.

  // Declare a map of prime factors raised to the power of the argument s.
  boost::unordered_map<std::uint32_t, local_complex_type> n_pow_s_prime_factor_map;

  // Generate a list of the first 300 prime numbers.
  static std::deque<std::uint32_t> prime_data;

  if(prime_data.empty())
  {
    ef::prime(static_cast<std::uint32_t>(300u), prime_data);
  }

  static const std::vector<std::uint32_t> primes(prime_data.begin(), prime_data.end());

  const local_real_type upper_limit { static_cast<std::uint32_t>(UINT32_C(1000100)) };

  using std::fabs;

  if(fabs(imag(s)) > upper_limit)
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
  constexpr double nd { static_cast<double>(std::numeric_limits<local_real_type>::digits10) * static_cast<double>(1.45) };

  using std::fabs;

  const double ni = static_cast<double>(static_cast<double>(1.10) * fabs(ef::to_double(imag(s))));

  const std::int32_t N { static_cast<std::int32_t>(static_cast<std::int64_t>(static_cast<double>(nd + ni))) };

  auto neg_term = ((N % static_cast<std::int32_t>(2)) == static_cast<std::int32_t>(0));

  const unsigned int
    two_n_minus_one
    {
      static_cast<unsigned int>
      (
        static_cast<std::int32_t>
        (
          static_cast<std::int32_t>(N + N) - static_cast<std::int32_t>(INT8_C(1))
        )
      )
    };

  local_real_type n_plus_j_minus_one_fact = boost::math::factorial<local_real_type>(two_n_minus_one);
  local_real_type four_pow_j              = Util::pown(local_real_type(static_cast<unsigned>(UINT8_C(4))), static_cast<std::int64_t>(N));
  local_real_type n_minus_j_fact          = ef::one<local_real_type>();
  local_real_type two_j_fact              = n_plus_j_minus_one_fact * static_cast<std::int32_t>(static_cast<std::int32_t>(INT8_C(2)) * N);

  local_real_type dn = (n_plus_j_minus_one_fact * four_pow_j) / (n_minus_j_fact * two_j_fact);

  local_complex_type jps = Util::j_pow_x(static_cast<std::uint32_t>(N), s, n_pow_s_prime_factor_map);

  local_complex_type zs = ((!neg_term) ? dn : -dn) / jps;

  for(auto   j  = static_cast<std::int32_t>(N - static_cast<std::int32_t>(INT8_C(1)));
              j >= static_cast<std::int32_t>(INT8_C(0));
            --j)
  {
    const auto j_is_zero = (j == static_cast<std::int32_t>(INT8_C(0)));

    const auto two_jp1_two_j =
      static_cast<std::int32_t>
      (
          static_cast<std::int32_t>((static_cast<std::int32_t>(INT8_C(2)) * j) + static_cast<std::int32_t>(INT8_C(1)))
        * static_cast<std::int32_t> (static_cast<std::int32_t>(INT8_C(2)) * (!j_is_zero ? j : static_cast<std::int32_t>(INT8_C(1))))
      );

    n_plus_j_minus_one_fact /= static_cast<std::int32_t>(N + j);
    four_pow_j              /= static_cast<std::int32_t>(INT8_C(4));
    n_minus_j_fact          *= static_cast<std::int32_t>(N - j);
    two_j_fact              /= two_jp1_two_j;

    dn += ((n_plus_j_minus_one_fact * four_pow_j) / (n_minus_j_fact * two_j_fact));

    if(!j_is_zero)
    {
      // Increment the zeta function sum.
      jps = Util::j_pow_x(j, s, n_pow_s_prime_factor_map);

      neg_term = (!neg_term);

      zs += ((!neg_term) ? dn : -dn) / jps;
    }
  }

  const local_complex_type one_minus_s { ef::one<local_real_type>() - s };

  static const local_complex_type two_cpx { local_real_type { static_cast<int>(INT8_C(2)) } };

  using std::pow;

  return zs / (dn * (ef::one<local_real_type>() - pow(two_cpx, one_minus_s)));
}

} // namespace detail

} // namespace zeta_detail

template<typename ComplexType>
auto riemann_zeta(const ComplexType& s) -> ComplexType
{
  return zeta_detail::detail::ZetaTemplate<ComplexType>(s);
}

#endif // ZETA_DETAIL_2024_03_17_H