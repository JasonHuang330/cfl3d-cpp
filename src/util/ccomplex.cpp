// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ccomplex.h"
#include <algorithm>
#include <cmath>

namespace ccomplex_ns {

// Purpose: To return a value of a or b, depending on whether the
// logical expression c is true or false
double ccvmgt(double& a, double& b, bool& c) {
    if (c) {
        return a;
    } else {
        return b;
    }
}

// Purpose: complex/real max
double ccmax(double& a1, double& a2) {
    return std::max(a1, a2);
}

// Purpose: complex/real max with 4 complex arguments
double ccmax4(double& a1, double& a2, double& a3, double& a4) {
    return std::max({a1, a2, a3, a4});
}

// Purpose: complex/real max with 1st argument complex, and 2nd argument real
double ccmaxcr(double& a1, double& a2) {
    return std::max(a1, a2);
}

// Purpose: complex/real max with 1st argument real, and 2nd argument complex
double ccmaxrc(double& a1, double& a2) {
    return std::max(a1, a2);
}

// Purpose: complex/real min
double ccmin(double& a1, double& a2) {
    return std::min(a1, a2);
}

// Purpose: complex/real min with 4 complex arguments
double ccmin4(double& a1, double& a2, double& a3, double& a4) {
    return std::min({a1, a2, a3, a4});
}

// Purpose: complex/real min with 8 complex arguments
double ccmin8(double& a1, double& a2, double& a3, double& a4,
              double& a5, double& a6, double& a7, double& a8) {
    return std::min({a1, a2, a3, a4, a5, a6, a7, a8});
}

// Purpose: complex/real min with 1st argument complex, and 2nd argument real
double ccmincr(double& a1, double& a2) {
    return std::min(a1, a2);
}

// Purpose: complex/real min with 1st argument real, and 2nd argument complex
double ccminrc(double& a1, double& a2) {
    return std::min(a1, a2);
}

// Purpose: complex/real absolute value
double ccabs(double& a) {
    return std::abs(a);
}

// Purpose: complex/real sign
// Fortran SIGN(a,b) = abs(a) * sign(b)
double ccsign(double& a, double& b) {
    return std::copysign(std::abs(a), b);
}

// Purpose: complex/real sign with 1st argument real, and 2nd argument complex
double ccsignrc(double& a, double& b) {
    return std::copysign(std::abs(a), b);
}

// Purpose: complex/real hyperbolic tangent
double cctanh(double& a) {
    return std::tanh(a);
}

// Purpose: complex arc tangent in correct quadrant
double ccatan2(double& csn, double& ccs) {
    return std::atan2(csn, ccs);
}

// Purpose: complex arc tangent (1st and 4th quadrants)
double ccatan(double& z) {
    return std::atan(z);
}

// Purpose: complex arccosine
double ccacos(double& z) {
    return std::acos(z);
}

// Purpose: complex arcsine
double ccasin(double& z) {
    return std::asin(z);
}

// Purpose: complex tangent
double cctan(double& z) {
    return std::tan(z);
}

// Purpose: complex log base 10
double cclog10(double& z) {
    return std::log10(z);
}

// Purpose: complex positive difference
// Fortran DIM(a1,a2) = max(a1-a2, 0)
double ccdim(double& a1, double& a2) {
    return std::max(a1 - a2, 0.0);
}

// Purpose: complex error function
double ccerf(double& a1) {
    return std::erf(a1);
}

} // namespace ccomplex_ns
