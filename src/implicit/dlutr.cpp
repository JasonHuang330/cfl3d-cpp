// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dlutr.h"

namespace dlutr_ns {

void dlutr(int& nvmax, int& n, int& nmax, int& il, int& iu,
           FortranArray2DRef<double> a,
           FortranArray2DRef<double> b,
           FortranArray2DRef<double> c)
{
    // Purpose: Perform the scalar tridiagonal (LU) decomposition.
    //
    // inversion of scalar tridiagonal...a,b,c are scalars
    // f is forcing function and solution is output in f
    // solution is by upper triangularization with unity diagonal
    // block inversions use nonpivoted lu decomposition
    // il and iu are starting and finishing indices
    // b and c are overloaded

    int il1 = il + 1;
    int i   = il;

    // l-u decomposition
    for (int izz = 1; izz <= n; izz++) {
        b(izz, i) = 1.e0 / b(izz, i);
    }

    if (i == iu) goto label_1030;

    // c=ainv*c
    for (int izz = 1; izz <= n; izz++) {
        c(izz, i) = b(izz, i) * c(izz, i);
    }

label_1030:
    // forward sweep
    for (i = il1; i <= iu; i++) {
        int ir = i - 1;
        // int it = i + 1;  // declared but not used in original

        // first row reduction
        for (int izz = 1; izz <= n; izz++) {
            b(izz, i) = b(izz, i) - a(izz, i) * c(izz, ir);

            // l-u decomposition
            b(izz, i) = 1.e0 / b(izz, i);
        }

        if (i == iu) continue;

        // c=ainv*c
        for (int izz = 1; izz <= n; izz++) {
            c(izz, i) = b(izz, i) * c(izz, i);
        }
    }

    return;
}

} // namespace dlutr_ns
