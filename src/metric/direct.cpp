// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "direct.h"
#include "ccomplex.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <cstring>

namespace direct_ns {

void direct(double& x5, double& x6, double& x7, double& x8,
            double& y5, double& y6, double& y7, double& y8,
            double& z5, double& z6, double& z7, double& z8,
            double& a1, double& a2, double& a3, int& imaxa,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim)
{
    double xxie = x8 - x7;
    double yxie = y8 - y7;
    double zxie = z8 - z7;
    double xeta = x6 - x5;
    double yeta = y6 - y5;
    double zeta = z6 - z5;
    a1 = yxie * zeta - zxie * yeta;
    a2 = zxie * xeta - xxie * zeta;
    a3 = xxie * yeta - yxie * xeta;
    double d = std::sqrt(a1 * a1 + a2 * a2 + a3 * a3);
    if (d <= 0.) {
        nou(4) = std::min(nou(4) + 1, ibufdim);
        std::snprintf(bou(nou(4), 4), 120,
            " WARNINING: a cell with zero area has been detected in subroutine direct");
        nou(4) = std::min(nou(4) + 1, ibufdim);
        std::snprintf(bou(nou(4), 4), 120, "  - severe problem");
        d = 1.;
    }
    a1 = a1 / d;
    a2 = a2 / d;
    a3 = a3 / d;
//
//   find coordinate direction with maximum area component
//
    double amax = ccomplex_ns::ccabs(a1);
    imaxa = 1;
    if (std::abs(a2) > amax) {
        amax  = ccomplex_ns::ccabs(a2);
        imaxa = 2;
    }
    if (std::abs(a3) > amax) {
        amax  = ccomplex_ns::ccabs(a3);
        imaxa = 3;
    }
    return;
}

} // namespace direct_ns
