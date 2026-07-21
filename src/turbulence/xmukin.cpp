// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "xmukin.h"
#include <cmath>

namespace xmukin_ns {

void xmukin(int& n, FortranArray1DRef<double> temp, FortranArray1DRef<double> visc, double& tinf)
{
//***********************************************************************
//     Purpose:  Computes Sutherland's formula.  Note that this routine
//     is only called by ctime.  (Mu is generally computed in-line
//     throughout the code using this formula.)
//     Also note that an older version of this routine used a linear
//     law for low temperatures, but it was not consistently used
//     throughout the whole code, and so has been taken out.
//***********************************************************************

    // common /fluid2/ pr, prt, cbar
    double cbar = (double)cmn_fluid2.cbar;

    // molecular viscosity
    double c2b  = cbar / tinf;
    double c2bp = c2b + 1.e0;

    for (int j = 1; j <= n; j++) {
        double t5   = temp(j);
        double t6   = std::sqrt(t5);
        visc(j)     = c2bp * t5 * t6 / (c2b + t5);
    }
}

} // namespace xmukin_ns
