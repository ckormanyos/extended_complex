///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 1999 - 2025.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef UTIL_2024_08_19_H
  #define UTIL_2024_08_19_H

  namespace util
  {
    template<typename NumericType>
    auto is_close_fraction(const NumericType& a,
                           const NumericType& b,
                           const NumericType& tol) noexcept -> bool
    {
      using std::fabs;
      using std::fpclassify;

      auto result_is_ok = bool { };

      if(fpclassify(b) == FP_ZERO)
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
  } // namespace util

#endif // UTIL_2024_08_19_H
