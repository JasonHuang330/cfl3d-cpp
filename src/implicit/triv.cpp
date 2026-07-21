// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "triv.h"

namespace triv_ns {

void triv(int& jdim, int& kdim, int& jl, int& ju, int& kl, int& ku,
          FortranArray2DRef<double> x, FortranArray2DRef<double> a,
          FortranArray2DRef<double> b, FortranArray2DRef<double> c,
          FortranArray2DRef<double> f)
{
    // Solve a scalar tridiagonal system of equations

    double z;
    int klp1, kupkl, i, i1, j;

    // Initialize first level
    for (j = jl; j <= ju; j++) {
        x(j, kl) = c(j, kl) / b(j, kl);
        f(j, kl) = f(j, kl) / b(j, kl);
    }

    // Forward sweep
    klp1 = kl + 1;
    for (i = klp1; i <= ku; i++) {
        for (j = jl; j <= ju; j++) {
            z = 1.0 / (b(j, i) - a(j, i) * x(j, i - 1));
            x(j, i) = c(j, i) * z;
            f(j, i) = (f(j, i) - a(j, i) * f(j, i - 1)) * z;
        }
    }

    // Back substitution
    kupkl = ku + kl;
    for (i1 = klp1; i1 <= ku; i1++) {
        i = kupkl - i1;
        for (j = jl; j <= ju; j++) {
            f(j, i) = f(j, i) - x(j, i) * f(j, i + 1);
        }
    }
}

} // namespace triv_ns
