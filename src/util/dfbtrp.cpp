// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dfbtrp.h"
#include "dfbtr.h"

namespace dfbtrp_ns {

void dfbtr(int& nvmax, int& n, int& nmax, int& il, int& iu,
           FortranArray2DRef<double> a, FortranArray2DRef<double> b,
           FortranArray2DRef<double> c, FortranArray2DRef<double> f)
{
    dfbtr_ns::dfbtr(nvmax, n, nmax, il, iu, a, b, c, f);
}

void dfbtrp(int& nvmax, int& n, int& nmax, int& il, int& iu,
            FortranArray2DRef<double> a, FortranArray2DRef<double> b,
            FortranArray2DRef<double> c, FortranArray2DRef<double> f,
            FortranArray2DRef<double> g, FortranArray2DRef<double> h)
{
//
//     Purpose:  Perform the back substitution for a periodic scalar
//     tridiagonal system of equations.
//
//      periodic solver  5-3-85
//
//      inversion of block tridiagonal...a,b,c are scalars
//      f is forcing function and solution is output in f
//      solution is by upper triangularization with unity diagonal
//      block inversions use nonpivoted lu decomposition
//      il and iu are starting and finishing indices
//      b,c,and e are overloaded
//
    int il1, is, iux, ir, it, iu2, il11, ii, ix, izz;

    il1 = il + 1;
    is  = il;

//      f=binv*f
    for (izz = 1; izz <= n; izz++) {
        f(izz, is) = b(izz, is) * (f(izz, is));
    }

//      forward sweep
    iux = iu;
    iux = iu - 1;
    for (is = il1; is <= iux; is++) {
        ir = is - 1;
        it = is + 1;
//      first row reduction
        for (izz = 1; izz <= n; izz++) {
            f(izz, is) = b(izz, is) * (f(izz, is) - a(izz, is) * f(izz, ir));
//
//      f=binv*f
//
        }
    }
    is = iu;
    ir = is - 1;
    it = is + 1;
//      first row reduction
    for (izz = 1; izz <= n; izz++) {
        f(izz, is) = f(izz, is) - a(izz, is) * f(izz, ir);
    }

    iu2 = iu - 2;
    for (ix = il; ix <= iu2; ix++) {
        for (izz = 1; izz <= n; izz++) {
            f(izz, is) = f(izz, is) - h(izz, ix) * f(izz, ix);
        }
    }

//      f=binv*f
    for (izz = 1; izz <= n; izz++) {
        f(izz, is) = b(izz, is) * (f(izz, is));
    }

//      back substitution
    iux = iu;
    iux = il1;
    for (ii = il1; ii <= iux; ii++) {
        is = il + iu - ii;
        it = is + 1;
        for (izz = 1; izz <= n; izz++) {
            f(izz, is) = f(izz, is) - c(izz, is) * f(izz, it);
        }
    }
    il11 = il1 + 1;
    for (ii = il11; ii <= iu; ii++) {
        is = il + iu - ii;
        it = is + 1;
        for (izz = 1; izz <= n; izz++) {
            f(izz, is) = f(izz, is) - c(izz, is) * f(izz, it) - g(izz, is) * f(izz, iu);
        }
    }
    return;
}

} // namespace dfbtrp_ns
