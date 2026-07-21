// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dfbtr.h"

namespace dfbtr_ns {

void dfbtr(int& nvmax, int& n, int& nmax, int& il, int& iu,
           FortranArray2DRef<double> a, FortranArray2DRef<double> b,
           FortranArray2DRef<double> c, FortranArray2DRef<double> f)
{
    int il1, is, ir, it, ii, izz;

    il1 = il + 1;
    is  = il;

    //     f=binv*f
    for (izz = 1; izz <= n; izz++) {
        f(izz, is) = b(izz, is) * f(izz, is);
    }

    //     forward sweep
    for (is = il1; is <= iu; is++) {
        ir = is - 1;
        it = is + 1;
        //     first row reduction
        //     f=binv*f
        for (izz = 1; izz <= n; izz++) {
            f(izz, is) = b(izz, is) * (f(izz, is) - a(izz, is) * f(izz, ir));
        }
    }

    //     back substitution
    for (ii = il1; ii <= iu; ii++) {
        is = il + iu - ii;
        it = is + 1;
        for (izz = 1; izz <= n; izz++) {
            f(izz, is) = f(izz, is) - c(izz, is) * f(izz, it);
        }
    }

    return;
}

} // namespace dfbtr_ns
