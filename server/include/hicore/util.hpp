#ifndef UTIL_HPP_
#define UTIL_HPP_
#include <cmath>
#include <limits>

static const double EPSILON = std::numeric_limits<double>::epsilon();

inline bool Equal(double a, double b) {
    return abs(a - b) <= ((abs(a) > abs(b) ? abs(b) : abs(a)) * EPSILON);
}

inline double Pow2(double x) { return x * x; }

#endif