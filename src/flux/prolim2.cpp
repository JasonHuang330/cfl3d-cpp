// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "prolim2.h"
#include "ccomplex.h"
#include "prolim.h"

namespace prolim2_ns {

void prolim(int& n, FortranArray1DRef<double> x1, FortranArray1DRef<double> x2, FortranArray1DRef<double> xc, int& leq)
{
    prolim_ns::prolim(n, x1, x2, xc, leq);
}

void prolim2(int& n, FortranArray1DRef<double> x1, FortranArray1DRef<double> x2, FortranArray1DRef<double> xc, int& leq)
{
    // data pmin,rhomin /0.01,0.01/
    static double pmin   = 0.01;
    static double rhomin = 0.01;

    if (leq == 1) {
        for (int izz = 1; izz <= n; izz++) {
            xc(izz) = ccomplex_ns::ccmax(xc(izz), rhomin);
            if ((float)xc(izz) == (float)rhomin) {
                x1(izz) = 0.;
                x2(izz) = 0.;
            }
        }
    } else if (leq == 5) {
        for (int izz = 1; izz <= n; izz++) {
            xc(izz) = ccomplex_ns::ccmax(xc(izz), pmin);
            if ((float)xc(izz) == (float)pmin) {
                x1(izz) = 0.;
                x2(izz) = 0.;
            }
        }
    }
}

} // namespace prolim2_ns
