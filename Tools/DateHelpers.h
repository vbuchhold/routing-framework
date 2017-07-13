#pragma once

#include <cstring>
#include <stdexcept>
#include <string>

#include "Tools/LexicalCast.h"
#include "Tools/StringHelpers.h"

// A day-of-week, such as 'Tuesday'.
enum class DayOfWeek {
  MONDAY    = 1,
  TUESDAY   = 2,
  WEDNESDAY = 3,
  THURSDAY  = 4,
  FRIDAY    = 5,
  SATURDAY  = 6,
  SUNDAY    = 7,
};

// Obtains the seconds since midnight from a text string such as 10:15:30.
template <typename StringT>
inline int parseTime(StringT& text) {
  const int len = str::length(text);
  if ((len == 5 || (len == 8 && text[5] == ':')) && text[2] == ':') {
    text[2] = '\0';
    text[5] = '\0';
    const int hrs = lexicalCast<int>(str::cStr(text));
    const int min = lexicalCast<int>(str::cStr(text) + 3);
    const int sec = len == 8 ? lexicalCast<int>(str::cStr(text) + 6) : 0;
    if (0 <= hrs && hrs < 24 && 0 <= min && min < 60 && 0 <= sec && sec < 60) {
      return hrs * 60 * 60 + min * 60 + sec;
    } else {
      text[2] = ':';
      if (len == 8)
        text[5] = ':';
    }
  }
  throw std::invalid_argument("text string not in format HH:mm:ss -- '" + std::string(text) + "'");
}

// Obtains an instance of DayOfWeek from a text string such as 'Tue'.
template <typename StringT>
inline DayOfWeek parseDayOfWeek(const StringT& text) {
  const char* const dayOfWeek = str::cStr(text);
  if (!std::strcmp(dayOfWeek, "Monday") || !std::strcmp(dayOfWeek, "Mon") ||
      !std::strcmp(dayOfWeek, "Montag") || !std::strcmp(dayOfWeek, "Mo"))
    return DayOfWeek::MONDAY;
  else if (!std::strcmp(dayOfWeek, "Tuesday") || !std::strcmp(dayOfWeek, "Tue") ||
      !std::strcmp(dayOfWeek, "Dienstag") || !std::strcmp(dayOfWeek, "Di"))
    return DayOfWeek::TUESDAY;
  else if (!std::strcmp(dayOfWeek, "Wednesday") || !std::strcmp(dayOfWeek, "Wed") ||
      !std::strcmp(dayOfWeek, "Mittwoch") || !std::strcmp(dayOfWeek, "Mi"))
    return DayOfWeek::WEDNESDAY;
  else if (!std::strcmp(dayOfWeek, "Thursday") || !std::strcmp(dayOfWeek, "Thu") ||
      !std::strcmp(dayOfWeek, "Donnerstag") || !std::strcmp(dayOfWeek, "Do"))
    return DayOfWeek::THURSDAY;
  else if (!std::strcmp(dayOfWeek, "Friday") || !std::strcmp(dayOfWeek, "Fri") ||
      !std::strcmp(dayOfWeek, "Freitag") || !std::strcmp(dayOfWeek, "Fr"))
    return DayOfWeek::FRIDAY;
  else if (!std::strcmp(dayOfWeek, "Saturday") || !std::strcmp(dayOfWeek, "Sat") ||
      !std::strcmp(dayOfWeek, "Samstag") || !std::strcmp(dayOfWeek, "Sa"))
    return DayOfWeek::SATURDAY;
  else if (!std::strcmp(dayOfWeek, "Sunday") || !std::strcmp(dayOfWeek, "Sun") ||
      !std::strcmp(dayOfWeek, "Sonntag") || !std::strcmp(dayOfWeek, "So"))
    return DayOfWeek::SUNDAY;
  throw std::invalid_argument("invalid day of week -- '" + std::string(text) + "'");
}

// Obtains the seconds since Monday midnight from two text strings such as 'Tue' and 10:15:30.
template <typename StringT1, typename StringT2>
inline int secondsSinceMonMidnight(const StringT1& dayOfWeek, StringT2& time) {
  return (static_cast<int>(parseDayOfWeek(dayOfWeek)) - 1) * 24 * 60 * 60 + parseTime(time);
}
