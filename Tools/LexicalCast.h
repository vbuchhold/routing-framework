#pragma once

#include <stdexcept>
#include <string>
#include <type_traits>

#include <csv.h>

#include "Tools/StringHelpers.h"

// Converts the specified string to type T.
template <typename T, typename StringT, typename = std::enable_if_t<!std::is_enum<T>::value>>
inline T lexicalCast(const StringT& string) {
  try {
    T val;
    io::detail::parse<io::throw_on_overflow>(const_cast<char*>(str::cStr(string)), val);
    return val;
  } catch (io::error::invalid_single_character& e) {
    const auto what = "'" + std::string(string) + "' is not a single character";
    throw std::invalid_argument(what);
  } catch (io::error::no_digit& e) {
    const auto what = "'" + std::string(string) + "' cannot be converted to an arithmetic type";
    throw std::invalid_argument(what);
  } catch (io::error::integer_overflow& e) {
    const auto what = "'" + std::string(string) + "' is outside the range of representable values";
    throw std::out_of_range(what);
  } catch (io::error::integer_underflow& e) {
    const auto what = "'" + std::string(string) + "' is outside the range of representable values";
    throw std::out_of_range(what);
  }
}

// Converts the specified string to the specified enum type.
template <typename EnumT, typename StringT>
inline std::enable_if_t<std::is_enum<EnumT>::value, EnumT> lexicalCast(const StringT& string) {
  return static_cast<EnumT>(lexicalCast<std::underlying_type_t<EnumT>>(string));
}
