// Basic types shared by the pyfda C++ core library (no Qt dependency).
#pragma once

#include <array>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace pyfda {

using cplx = std::complex<double>;
using Vec = std::vector<double>;
using CVec = std::vector<cplx>;

/// Second-order sections: each row is [b0, b1, b2, a0, a1, a2] (like scipy)
using Sos = std::vector<std::array<double, 6>>;

/// Zeros, poles and gain
struct Zpk {
    CVec z;
    CVec p;
    double k = 1.0;
};

/// Transfer function coefficients (numerator b, denominator a)
struct Ba {
    Vec b;
    Vec a;
};

/// Exception thrown for invalid design parameters, messages are shown to the user
class DesignError : public std::runtime_error {
public:
    explicit DesignError(const std::string &msg) : std::runtime_error(msg) {}
};

constexpr double PI = 3.14159265358979323846;

}  // namespace pyfda
