#pragma once

// Returns its argument, doing nothing with it at all. This function helps to avoid odr-using
// constant static data members. Constant static data members that are not odr-used in the program
// aren't required to be defined in a namespace scope, which is a good thing. In C++17, the better
// solution is to use inline variables. However, this is C++14 code.
template <typename T>
T use(T arg) { return arg; }

// Marks its arguments as possibly unused to suppress compiler warnings.
template <typename ...Types>
void unused(const Types&...) {}
