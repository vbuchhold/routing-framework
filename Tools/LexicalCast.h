#pragma once

#include <stdexcept>
#include <string>
#include <type_traits>

// Converts the specified string to type T.
template <typename T, typename = std::enable_if_t<!std::is_enum<T>::value>>
inline T lexicalCast(const std::string& str);

// Converts the specified string to the specified integral type, using the specified base.
template <typename IntegralT>
inline IntegralT lexicalCast(const std::string& str, const int base);

// Converts the specified string to type int, using the specified base.
template <>
inline int lexicalCast(const std::string& str, const int base) {
  size_t firstUnconverted;
  int res = std::stoi(str, &firstUnconverted, base);
  if (firstUnconverted != str.size())
    throw std::invalid_argument("'" + str + "' cannot be converted to type int");
  return res;
}

// Converts the specified string to type char.
template <>
inline char lexicalCast(const std::string& str) {
  if (str.size() != 1)
    throw std::invalid_argument("'" + str + "' cannot be converted to type char");
  return str[0];
}

// Converts the specified string to type int.
template <>
inline int lexicalCast(const std::string& str) {
  return lexicalCast<int>(str, 10);
}

// Converts the specified string to type float.
template <>
inline float lexicalCast(const std::string& str) {
  size_t firstUnconverted;
  float res = std::stof(str, &firstUnconverted);
  if (firstUnconverted != str.size())
    throw std::invalid_argument("'" + str + "' cannot be converted to type float");
  return res;
}

// Converts the specified string to type double.
template <>
inline double lexicalCast(const std::string& str) {
  size_t firstUnconverted;
  double res = std::stod(str, &firstUnconverted);
  if (firstUnconverted != str.size())
    throw std::invalid_argument("'" + str + "' cannot be converted to type double");
  return res;
}

// Converts the specified string to type long double.
template <>
inline long double lexicalCast(const std::string& str) {
  size_t firstUnconverted;
  long double res = std::stold(str, &firstUnconverted);
  if (firstUnconverted != str.size())
    throw std::invalid_argument("'" + str + "' cannot be converted to type long double");
  return res;
}

// Converts the specified string to type std::string.
template <>
inline std::string lexicalCast(const std::string& str) {
  return str;
}

// Converts the specified string to the specified enum type.
template <typename EnumT>
inline std::enable_if_t<std::is_enum<EnumT>::value, EnumT> lexicalCast(const std::string& str) {
  return static_cast<EnumT>(lexicalCast<int>(str));
}
