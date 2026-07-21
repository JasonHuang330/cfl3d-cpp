// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dlutrp.h"
#include "dlutr.h"

namespace dlutrp_ns {

void dlutr(int& nvmax, int& n, int& nmax, int& il, int& iu,
           FortranArray2DRef<double> a, FortranArray2DRef<double> b,
           FortranArray2DRef<double> c)
{
    dlutr_ns::dlutr(nvmax, n, nmax, il, iu, a, b, c);
}

void dlutrp(int& nvmax, int& n, int& nmax, int& il, int& iu,
            FortranArray2DRef<double> a, FortranArray2DRef<double> b,
            FortranArray2DRef<double> c, FortranArray2DRef<double> g,
            FortranArray2DRef<double> h)
{
    //
    //     periodic solver  5-3-85
    //
    //     inversion of block tridiagonal...a,b,c are scalars
    //     f is forcing function and solution is output in f
    //     solution is by upper triangularization with unity diagonal
    //     block inversions use nonpivoted lu decomposition
    //     il and iu are starting and finishing indices
    //     b,c,and e are overloaded
    //
    //     set g and h
    //
    for (int izz = 1; izz <= n; izz++) {
        g(izz, il) = a(izz, il);
        h(izz, il) = c(izz, iu);
    }
    //
    int il1 = il + 1;
    int i   = il;
    //
    //     l-u decomposition
    //
    for (int izz = 1; izz <= n; izz++) {
        b(izz, i) = 1.e0 / b(izz, i);
    }
    //
    if (i == iu) goto label_1030;
    //
    //     c=binv*c
    //     g=binv*g
    //
    for (int izz = 1; izz <= n; izz++) {
        c(izz, i) = b(izz, i) * (c(izz, i));
        g(izz, i) = b(izz, i) * (g(izz, i));
    }
label_1030:
    //     forward sweep
    {
        int iu2 = iu - 2;
        if (iu2 == il) goto label_101;
        for (i = il1; i <= iu2; i++) {
            int ir = i - 1;
            int it = i + 1;
            (void)it;
            //     row reduction
            for (int izz = 1; izz <= n; izz++) {
                b(izz, i)  =  b(izz, i)  - a(izz, i) * c(izz, ir);
                g(izz, i)  = -a(izz, i)  * g(izz, ir);
                h(izz, i)  = -h(izz, ir) * c(izz, ir);
                b(izz, iu) =  b(izz, iu) - h(izz, ir) * g(izz, ir);
                //
                //     l-u decomposition
                //
                b(izz, i) = 1.e0 / b(izz, i);
                //
                //     c=binv*c
                //
                c(izz, i) = b(izz, i) * (c(izz, i));
                g(izz, i) = b(izz, i) * (g(izz, i));
            }
        }
    }
label_101:
    //
    i      = iu - 1;
    {
        int ir = i - 1;
        int it = i + 1;
        (void)it;
        //     row reduction
        for (int izz = 1; izz <= n; izz++) {
            b(izz, i)  = b(izz, i)  - a(izz, i) * c(izz, ir);
            c(izz, i)  = c(izz, i)  - a(izz, i) * g(izz, ir);
            //     ************************** 9-29-86
            a(izz, iu) = a(izz, iu) - h(izz, ir) * c(izz, ir);
            b(izz, iu) = b(izz, iu) - h(izz, ir) * g(izz, ir);
            //
            //     l-u decomposition
            //
            b(izz, i) = 1.e0 / b(izz, i);
            //
            //     c=binv*c
            //
            c(izz, i) = b(izz, i) * (c(izz, i));
        }
    }
    i      = iu;
    {
        int ir = i - 1;
        int it = i + 1;
        (void)it;
        //     row reduction
        for (int izz = 1; izz <= n; izz++) {
            b(izz, i) = b(izz, i) - a(izz, i) * c(izz, ir);
            //
            //     l-u decomposition
            //
            b(izz, i) = 1.e0 / b(izz, i);
        }
    }
    return;
}

} // namespace dlutrp_ns
