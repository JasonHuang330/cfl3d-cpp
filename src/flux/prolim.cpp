// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "prolim.h"
#include "ccomplex.h"
#include "runtime/common_blocks.h"

namespace prolim_ns {

void prolim(int& n, FortranArray1DRef<double> x1, FortranArray1DRef<double> x2, FortranArray1DRef<double> xc, int& leq)
{
//***********************************************************************
//     Purpose: to limit static pressure and static density to be less
//              than or equal to the stagnation value, and greater than
//              or equal to a minimum value set via the data statement.
//              For points that have the pressure/density so limited,
//              the difference stencil is also make locally first order
//              accurate.
//***********************************************************************

    // COMMON /fluid/
    // (not directly used in this routine but declared in Fortran)

    // COMMON /ivals/
    float& pt0   = cmn_ivals.pt0;
    float& rhot0 = cmn_ivals.rhot0;

    // data pmin,rhomin /0.1,0.1/
    static const double pmin   = 0.1;
    static const double rhomin = 0.1;

    if (leq == 1) {
        for (int izz = 1; izz <= n; izz++) {
            double rhomin_d = rhomin;
            double rhot0_d  = (double)rhot0;
            xc(izz) = ccomplex_ns::ccmax(xc(izz), rhomin_d);
            xc(izz) = ccomplex_ns::ccmin(xc(izz), rhot0_d);
            if ((float)(xc(izz) - rhomin) * (float)(xc(izz) - rhot0) == 0.f) {
                x1(izz) = 0.;
                x2(izz) = 0.;
            }
        }
    } else if (leq == 5) {
        for (int izz = 1; izz <= n; izz++) {
            double pmin_d = pmin;
            double pt0_d  = (double)pt0;
            xc(izz) = ccomplex_ns::ccmax(xc(izz), pmin_d);
            xc(izz) = ccomplex_ns::ccmin(xc(izz), pt0_d);
            if ((float)(xc(izz) - pmin) * (float)(xc(izz) - pt0) == 0.f) {
                x1(izz) = 0.;
                x2(izz) = 0.;
            }
        }
    }
}

} // namespace prolim_ns
