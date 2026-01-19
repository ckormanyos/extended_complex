///////////////////////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 1999 - 2026.
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef UTIL_2024_08_19_H
  #define UTIL_2024_08_19_H

  #include <cmath>
  #include <sstream>

  namespace util
  {
    template<typename FloatType>
    auto my_lexical_cast(const char* p_str) -> FloatType
    {
      std::stringstream strm { };

      strm << p_str;

      using float_type = FloatType;

      float_type flt { };

      strm >> flt;

      return flt;
    };

    template<typename NumericType>
    auto is_close_fraction(const NumericType& a,
                           const NumericType& b,
                           const NumericType& tol) -> bool
    {
      using std::fabs;
      using std::fpclassify;

      using numeric_type = NumericType;

      const numeric_type delta { (fpclassify(b) == FP_ZERO) ? fabs(a - b) : fabs(1 - (a / b)) };

      const bool result_is_ok { (delta < tol) };

      return result_is_ok;
    }
  } // namespace util

#endif // UTIL_2024_08_19_H
