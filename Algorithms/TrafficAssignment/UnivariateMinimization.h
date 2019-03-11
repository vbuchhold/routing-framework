#pragma once

#include <cassert>

// An implementation of the bisection method, also known as Bolzano search. Returns the minimum of
// a ditonic function in the interval [a, b], with a tolerance of +/- epsilon. Note that the first
// parameter is the derivative of the function to be minimized.
template <typename DerivativeT>
inline double bisectionMethod(DerivativeT derivative, double a, double b, double epsilon = 1e-15) {
  assert(a <= b);
  assert(epsilon > 0);

  if (derivative(a) >= 0)
    return a;
  if (derivative(b) <= 0)
    return b;

  // Repeatedly halve the current interval until it is no greater than 2 * epsilon.
  while (b - a > 2 * epsilon) {
    const auto midpoint = (b + a) / 2;
    if (derivative(midpoint) <= 0)
      a = midpoint;
    else
      b = midpoint;
  }
  return (b + a) / 2;
}
