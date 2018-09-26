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
  // Constructs the color black.
  constexpr Color() : redValue(0), greenValue(0), blueValue(0), alphaValue(255) {}

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
constexpr Color KIT_GREEN    = Color(  0, 150, 130);
constexpr Color KIT_GREEN_70 = Color( 77, 182, 168);
constexpr Color KIT_GREEN_50 = Color(127, 202, 192);
constexpr Color KIT_GREEN_30 = Color(178, 223, 217);
constexpr Color KIT_GREEN_15 = Color(217, 239, 236);

constexpr Color KIT_BLUE    = Color( 70, 100, 170);
constexpr Color KIT_BLUE_70 = Color(126, 147, 196);
constexpr Color KIT_BLUE_50 = Color(162, 177, 212);
constexpr Color KIT_BLUE_30 = Color(199, 208, 229);
constexpr Color KIT_BLUE_15 = Color(227, 232, 242);

constexpr Color KIT_BLACK    = Color(  0,   0,   0);
constexpr Color KIT_BLACK_70 = Color( 77,  77,  77);
constexpr Color KIT_BLACK_50 = Color(127, 127, 127);
constexpr Color KIT_BLACK_30 = Color(178, 178, 178);
constexpr Color KIT_BLACK_15 = Color(217, 217, 217);

constexpr Color KIT_PALEGREEN    = Color(140, 182,  60);
constexpr Color KIT_PALEGREEN_70 = Color(174, 204, 118);
constexpr Color KIT_PALEGREEN_50 = Color(197, 218, 157);
constexpr Color KIT_PALEGREEN_30 = Color(220, 233, 196);
constexpr Color KIT_PALEGREEN_15 = Color(238, 244, 226);

constexpr Color KIT_YELLOW    = Color(252, 229,   0);
constexpr Color KIT_YELLOW_70 = Color(253, 237,  77);
constexpr Color KIT_YELLOW_50 = Color(253, 242, 127);
constexpr Color KIT_YELLOW_30 = Color(254, 247, 178);
constexpr Color KIT_YELLOW_15 = Color(155, 251, 217);

constexpr Color KIT_ORANGE    = Color(223, 155,  27);
constexpr Color KIT_ORANGE_70 = Color(233, 185,  95);
constexpr Color KIT_ORANGE_50 = Color(239, 205, 141);
constexpr Color KIT_ORANGE_30 = Color(245, 225, 186);
constexpr Color KIT_ORANGE_15 = Color(250, 240, 221);

constexpr Color KIT_BROWN    = Color(167, 130,  46);
constexpr Color KIT_BROWN_70 = Color(193, 167, 108);
constexpr Color KIT_BROWN_50 = Color(211, 192, 150);
constexpr Color KIT_BROWN_30 = Color(228, 217, 192);
constexpr Color KIT_BROWN_15 = Color(242, 236, 224);

constexpr Color KIT_RED    = Color(162,  34,  35);
constexpr Color KIT_RED_70 = Color(190, 100, 101);
constexpr Color KIT_RED_50 = Color(208, 144, 145);
constexpr Color KIT_RED_30 = Color(227, 188, 189);
constexpr Color KIT_RED_15 = Color(241, 222, 222);

constexpr Color KIT_LILAC    = Color(163,  16, 124);
constexpr Color KIT_LILAC_70 = Color(190,  87, 163);
constexpr Color KIT_LILAC_50 = Color(209, 135, 189);
constexpr Color KIT_LILAC_30 = Color(227, 183, 215);
constexpr Color KIT_LILAC_15 = Color(241, 219, 235);

constexpr Color KIT_CYANBLUE    = Color( 35, 161, 224);
constexpr Color KIT_CYANBLUE_70 = Color(101, 189, 233);
constexpr Color KIT_CYANBLUE_50 = Color(145, 208, 239);
constexpr Color KIT_CYANBLUE_30 = Color(189, 227, 246);
constexpr Color KIT_CYANBLUE_15 = Color(222, 241, 250);

// Some predefined color schemes.
constexpr std::array<Color, 10> KIT_SCHEME = {{
    KIT_GREEN, KIT_BLUE, KIT_BLACK, KIT_PALEGREEN, KIT_YELLOW,
    KIT_ORANGE, KIT_BROWN, KIT_RED, KIT_LILAC, KIT_CYANBLUE}};

constexpr std::array<Color, 10> KIT_SCHEME_70 = {{
    KIT_GREEN_70, KIT_BLUE_70, KIT_BLACK_70, KIT_PALEGREEN_70, KIT_YELLOW_70,
    KIT_ORANGE_70, KIT_BROWN_70, KIT_RED_70, KIT_LILAC_70, KIT_CYANBLUE_70}};

constexpr std::array<Color, 10> KIT_SCHEME_50 = {{
    KIT_GREEN_50, KIT_BLUE_50, KIT_BLACK_50, KIT_PALEGREEN_50, KIT_YELLOW_50,
    KIT_ORANGE_50, KIT_BROWN_50, KIT_RED_50, KIT_LILAC_50, KIT_CYANBLUE_50}};

constexpr std::array<Color, 10> KIT_SCHEME_30 = {{
    KIT_GREEN_30, KIT_BLUE_30, KIT_BLACK_30, KIT_PALEGREEN_30, KIT_YELLOW_30,
    KIT_ORANGE_30, KIT_BROWN_30, KIT_RED_30, KIT_LILAC_30, KIT_CYANBLUE_30}};

constexpr std::array<Color, 10> KIT_SCHEME_15 = {{
    KIT_GREEN_15, KIT_BLUE_15, KIT_BLACK_15, KIT_PALEGREEN_15, KIT_YELLOW_15,
    KIT_ORANGE_15, KIT_BROWN_15, KIT_RED_15, KIT_LILAC_15, KIT_CYANBLUE_15}};

constexpr std::array<Color, 9> REDS_9CLASS = {{
    Color(0xfff5f0), Color(0xfee0d2), Color(0xfcbba1),
    Color(0xfc9272), Color(0xfb6a4a), Color(0xef3b2c),
    Color(0xcb181d), Color(0xa50f15), Color(0x67000d)}};
