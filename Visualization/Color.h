#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>

// A color in the RGB color space. Every color has an implicit alpha value of 1.0 or an explicit
// one provided in the constructor. The color and alpha components can be specified by values in
// the range 0.0-1.0 or 0-255.
class Color {
 public:
  // Constructs a color with the specified red, green, blue and alpha values in the range 0-255.
  constexpr Color(const int red, const int green, const int blue, const int alpha = 255)
      : redValue(red), greenValue(green), blueValue(blue), alphaValue(alpha) {
    assert(red >= 0); assert(red <= 255);
    assert(green >= 0); assert(green <= 255);
    assert(blue >= 0); assert(blue <= 255);
    assert(alpha >= 0); assert(alpha <= 255);
  }

  // Constructs a color with the specified hexadecimal representation.
  explicit constexpr Color(const uint32_t rgba, const bool hasAlpha = false)
      : redValue(rgba >> 16),
        greenValue(rgba >> 8),
        blueValue(rgba),
        alphaValue(hasAlpha ? rgba >> 24 : 255) {}

  // Constructs a color with the specified red, green, blue and alpha values in the range 0.0-1.0.
  constexpr Color(const double red, const double green, const double blue, const double alpha = 1.0)
      : redValue(std::round(red * 255)),
        greenValue(std::round(green * 255)),
        blueValue(std::round(blue * 255)),
        alphaValue(std::round(alpha * 255)) {
    assert(red >= 0.0); assert(red <= 1.0);
    assert(green >= 0.0); assert(green <= 1.0);
    assert(blue >= 0.0); assert(blue <= 1.0);
    assert(alpha >= 0.0); assert(alpha <= 1.0);
  }

  // Returns the red value in the range 0-255.
  uint8_t red() const {
    return redValue;
  }

  // Returns the green value in the range 0-255.
  uint8_t green() const {
    return greenValue;
  }

  // Returns the blue value in the range 0-255.
  uint8_t blue() const {
    return blueValue;
  }

  // Returns the alpha value in the range 0-255.
  uint8_t alpha() const {
    return alphaValue;
  }

  // Returns the red value in the range 0.0-1.0.
  double normalizedRed() const {
    return redValue / 255.0;
  }

  // Returns the green value in the range 0.0-1.0.
  double normalizedGreen() const {
    return greenValue / 255.0;
  }

  // Returns the blue value in the range 0.0-1.0.
  double normalizedBlue() const {
    return blueValue / 255.0;
  }

  // Returns the alpha value in the range 0.0-1.0.
  double normalizedAlpha() const {
    return alphaValue / 255.0;
  }

 private:
  uint8_t redValue;
  uint8_t greenValue;
  uint8_t blueValue;
  uint8_t alphaValue;
};

// Some predefined colors.
constexpr Color KIT_GREEN = Color(0, 150, 130);
constexpr Color KIT_BLUE = Color(70, 100, 170);
constexpr Color KIT_BLACK = Color(0, 0, 0);
constexpr Color KIT_PALEGREEN = Color(130, 190, 60);
constexpr Color KIT_YELLOW = Color(250, 230, 20);
constexpr Color KIT_ORANGE = Color(220, 160, 30);
constexpr Color KIT_BROWN = Color(160, 130, 50);
constexpr Color KIT_RED = Color(160, 30, 40);
constexpr Color KIT_LILAC = Color(160, 0, 120);
constexpr Color KIT_CYANBLUE = Color(80, 170, 230);

// Some predefined color schemes.
constexpr std::array<Color, 9> REDS_9CLASS = {
  Color(0xfff5f0), Color(0xfee0d2), Color(0xfcbba1),
  Color(0xfc9272), Color(0xfb6a4a), Color(0xef3b2c),
  Color(0xcb181d), Color(0xa50f15), Color(0x67000d) };
