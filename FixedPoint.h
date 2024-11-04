// FixedPoint.h
#ifndef FIXEDPOINT_H
#define FIXEDPOINT_H

#include <ap_fixed.h>
#include <cmath>

// Define a 32-bit fixed-point type with 16 integer bits and 16 fractional bits
typedef ap_fixed<32, 16> fixed_t;

// Convert double to fixed-point
inline fixed_t double_to_fixed(double value) {
    return fixed_t(value);
}

// Convert fixed-point to double
inline double fixed_to_double(fixed_t value) {
    return value.to_double();
}

// Fixed-point multiplication
inline fixed_t fixed_mul(fixed_t a, fixed_t b) {
    return a * b;
}

// Fixed-point division
inline fixed_t fixed_div(fixed_t a, fixed_t b) {
    if (b == 0) {
        // Division by zero
        return 0;
    }
    return a / b;
}

// Fixed-point absolute value
inline fixed_t fixed_abs(fixed_t x) {
    return (x < 0) ? fixed_t(-x) : x;
}

#endif // FIXEDPOINT_H
