// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rie1de.h"
#include "rie1d.h"
#include "ccomplex.h"
#include <cmath>

namespace rie1de_ns {

void rie1d(int& jvdim, FortranArray2DRef<double> t, int& jv, double& cl)
{
    rie1d_ns::rie1d(jvdim, t, jv, cl);
}

void rie1de(int& jvdim, FortranArray2DRef<double> t, int& jv)
{
    // unsteady quasi-1-d characteristic boundary conditions
    // for inflow/outflow

    double gamma = (double)cmn_fluid.gamma;
    double gm1   = (double)cmn_fluid.gm1;

    // Loop 3117
    for (int izz = 1; izz <= jv; izz++) {
        t(izz, 10) = gamma * t(izz, 5) / t(izz, 1);
    }

    double x2gm1 = 2.e0 / gm1;

    // Loop 3118
    for (int izz = 1; izz <= jv; izz++) {
        t(izz, 11) = t(izz, 6) * t(izz, 2) + t(izz, 7) * t(izz, 3) + t(izz, 8) * t(izz, 4);
        t(izz, 12) = std::sqrt(t(izz, 10));
        t(izz, 12) = t(izz, 11) + x2gm1 * t(izz, 12);
        t(izz, 13) = t(izz, 6) * t(izz, 21) + t(izz, 7) * t(izz, 22)
                   + t(izz, 8) * t(izz, 23);
        t(izz, 14) = t(izz, 13) - x2gm1 * t(izz, 24);
    }

    double x4gm1 = 0.25e0 * gm1;

    // Loop 3119
    for (int izz = 1; izz <= jv; izz++) {
        t(izz, 15) = .5e0 * (t(izz, 12) + t(izz, 14));
        t(izz, 16) = x4gm1 * (t(izz, 12) - t(izz, 14));
        t(izz, 17) = t(izz, 15) - t(izz, 13);
        // put unsteady metrics into t(5)
        t(izz, 5)  = t(izz, 20);
        t(izz, 18) = t(izz, 21) + t(izz, 6) * t(izz, 17);
        t(izz, 19) = t(izz, 22) + t(izz, 7) * t(izz, 17);
        t(izz, 20) = t(izz, 23) + t(izz, 8) * t(izz, 17);
    }

    // Loop 3120
    for (int izz = 1; izz <= jv; izz++) {
        t(izz, 13) = t(izz, 15) + t(izz, 5);
        bool cond = ((float)t(izz, 13) >= 0.e0);
        double a17a = t(izz, 15) - t(izz, 11);
        double b17  = t(izz, 17);
        t(izz, 17) = ccomplex_ns::ccvmgt(a17a, b17, cond);
        double a18a = t(izz, 2) + t(izz, 6) * t(izz, 17);
        double b18  = t(izz, 18);
        t(izz, 18) = ccomplex_ns::ccvmgt(a18a, b18, cond);
        double a19a = t(izz, 3) + t(izz, 7) * t(izz, 17);
        double b19  = t(izz, 19);
        t(izz, 19) = ccomplex_ns::ccvmgt(a19a, b19, cond);
        double a20a = t(izz, 4) + t(izz, 8) * t(izz, 17);
        double b20  = t(izz, 20);
        t(izz, 20) = ccomplex_ns::ccvmgt(a20a, b20, cond);
    }

    double xgm1 = 1.e0 / gm1;

    // Loop 3121
    for (int izz = 1; izz <= jv; izz++) {
        double rmloc2 = (t(izz, 18) * t(izz, 18) + t(izz, 19) * t(izz, 19) + t(izz, 20) * t(izz, 20))
                      / (t(izz, 16) * t(izz, 16));
        t(izz, 5) = t(izz, 25) / (gamma * std::pow(1.0 + 0.5 * gm1 * rmloc2, gamma / gm1));
        t(izz, 1) = gamma * t(izz, 5) / (t(izz, 16) * t(izz, 16));
        t(izz, 2) = t(izz, 18);
        t(izz, 3) = t(izz, 19);
        t(izz, 4) = t(izz, 20);
    }
}

} // namespace rie1de_ns
