// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rotatp.h"
#include <cmath>

namespace rotatp_ns {

void rotatp(int& mdim, int& ndim, int& jmax1, int& kmax1, int& msub1, int& l,
            FortranArray3DRef<double> x1, FortranArray3DRef<double> y1,
            FortranArray3DRef<double> z1, FortranArray2DRef<double> dthetx,
            FortranArray2DRef<double> dthety, FortranArray2DRef<double> dthetz,
            FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig,
            FortranArray1DRef<double> zorig, int& mbl, int& nn, int& intmx,
            int& int_)
{
    double xorg, yorg, zorg;
    double dthx, dthy, dthz;
    double ca, sa;
    double xtemp, ytemp, ztemp;

    xorg = xorig(mbl);
    yorg = yorig(mbl);
    zorg = zorig(mbl);
    dthx = dthetx(int_, l);
    dthy = dthety(int_, l);
    dthz = dthetz(int_, l);

    if (std::abs((float)dthx) > 0.) {
        //
        //  rotation about a line parallel to the x-axis
        //
        ca = std::cos(dthx);
        sa = std::sin(dthx);
        for (int j = 1; j <= jmax1; j++) {
            for (int k = 1; k <= kmax1; k++) {
                ytemp = (y1(j,k,l) - yorg)*ca - (z1(j,k,l) - zorg)*sa + yorg;
                ztemp = (y1(j,k,l) - yorg)*sa + (z1(j,k,l) - zorg)*ca + zorg;
                y1(j,k,l) = ytemp;
                z1(j,k,l) = ztemp;
            }
        }
    }

    if (std::abs((float)dthy) > 0.) {
        //
        //  rotation about a line parallel to the y-axis
        //
        ca = std::cos(dthy);
        sa = std::sin(dthy);
        for (int j = 1; j <= jmax1; j++) {
            for (int k = 1; k <= kmax1; k++) {
                xtemp =  (x1(j,k,l) - xorg)*ca + (z1(j,k,l) - zorg)*sa + xorg;
                ztemp = -(x1(j,k,l) - xorg)*sa + (z1(j,k,l) - zorg)*ca + zorg;
                x1(j,k,l) = xtemp;
                z1(j,k,l) = ztemp;
            }
        }
    }

    if (std::abs((float)dthz) > 0.) {
        //
        //  rotation about a line parallel to the z-axis
        //
        ca = std::cos(dthz);
        sa = std::sin(dthx);  // NOTE: original Fortran uses dthx here (preserved faithfully)
        for (int j = 1; j <= jmax1; j++) {
            for (int k = 1; k <= kmax1; k++) {
                xtemp = (x1(j,k,l) - xorg)*ca - (y1(j,k,l) - yorg)*sa + xorg;
                ytemp = (x1(j,k,l) - xorg)*sa + (y1(j,k,l) - yorg)*ca + yorg;
                x1(j,k,l) = xtemp;
                y1(j,k,l) = ytemp;
            }
        }
    }

    return;
}

} // namespace rotatp_ns
