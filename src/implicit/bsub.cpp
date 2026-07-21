// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bsub.h"

namespace bsub_ns {

void bsub(int& id1, int& id2, FortranArray4DRef<double> a, FortranArray4DRef<double> b, FortranArray4DRef<double> c, FortranArray3DRef<double> f, int& i1, int& i2, int& il, int& iu)
{
    // Purpose: Performs the back substitution for a block 5x5 tridi-
    // agonal matrix equation solution. The vectorization is over points
    // i1-i2 and the tridiagonal matrix equation spans points il-iu.

    int il1, is, ir, it, iqq, m, i;
    double t1, t2, t3, t4;

    // solve matrix equation
    il1 = il + 1;
    is  = il;

    // f=binv*f
    for (i = i1; i <= i2; i++) {
        t1        = b(i,is,1,1)*(f(i,is,1));
        t2        = b(i,is,2,2)*(f(i,is,2)-b(i,is,2,1)*t1);
        t3        = b(i,is,3,3)*(f(i,is,3)-b(i,is,3,1)*t1-b(i,is,3,2)*t2);
        t4        = b(i,is,4,4)*(f(i,is,4)-b(i,is,4,1)*t1-b(i,is,4,2)*t2
                                            -b(i,is,4,3)*t3);
        f(i,is,5) = b(i,is,5,5)*(f(i,is,5)-b(i,is,5,1)*t1-b(i,is,5,2)*t2
                                            -b(i,is,5,3)*t3-b(i,is,5,4)*t4);
        f(i,is,4) = t4-b(i,is,4,5)*f(i,is,5);
        f(i,is,3) = t3-b(i,is,3,5)*f(i,is,5)-b(i,is,3,4)*f(i,is,4);
        f(i,is,2) = t2-b(i,is,2,5)*f(i,is,5)-b(i,is,2,4)*f(i,is,4)
                      -b(i,is,2,3)*f(i,is,3);
        f(i,is,1) = t1-b(i,is,1,5)*f(i,is,5)-b(i,is,1,4)*f(i,is,4)
                      -b(i,is,1,3)*f(i,is,3)-b(i,is,1,2)*f(i,is,2);
    }

    // forward sweep
    for (is = il1; is <= iu; is++) {
        ir = is - 1;
        it = is + 1;
        // first row reduction
        for (m = 1; m <= 5; m++) {
            for (i = i1; i <= i2; i++) {
                f(i,is,m) = f(i,is,m)-a(i,is,m,1)*f(i,ir,1)
                                      -a(i,is,m,2)*f(i,ir,2)
                                      -a(i,is,m,3)*f(i,ir,3)
                                      -a(i,is,m,4)*f(i,ir,4)
                                      -a(i,is,m,5)*f(i,ir,5);
            }
        }

        // f=binv*f
        for (i = i1; i <= i2; i++) {
            t1        = b(i,is,1,1)*(f(i,is,1));
            t2        = b(i,is,2,2)*(f(i,is,2)-b(i,is,2,1)*t1);
            t3        = b(i,is,3,3)*(f(i,is,3)-b(i,is,3,1)*t1-b(i,is,3,2)*t2);
            t4        = b(i,is,4,4)*(f(i,is,4)-b(i,is,4,1)*t1-b(i,is,4,2)*t2
                                                -b(i,is,4,3)*t3);
            f(i,is,5) = b(i,is,5,5)*(f(i,is,5)-b(i,is,5,1)*t1-b(i,is,5,2)*t2
                                                -b(i,is,5,3)*t3-b(i,is,5,4)*t4);
            f(i,is,4) = t4-b(i,is,4,5)*f(i,is,5);
            f(i,is,3) = t3-b(i,is,3,5)*f(i,is,5)-b(i,is,3,4)*f(i,is,4);
            f(i,is,2) = t2-b(i,is,2,5)*f(i,is,5)-b(i,is,2,4)*f(i,is,4)
                          -b(i,is,2,3)*f(i,is,3);
            f(i,is,1) = t1-b(i,is,1,5)*f(i,is,5)-b(i,is,1,4)*f(i,is,4)
                          -b(i,is,1,3)*f(i,is,3)-b(i,is,1,2)*f(i,is,2);
        }
    }

    // back substitution
    for (iqq = il1; iqq <= iu; iqq++) {
        is = il + iu - iqq;
        it = is + 1;
        for (m = 1; m <= 5; m++) {
            for (i = i1; i <= i2; i++) {
                f(i,is,m) =  f(i,is,m)-c(i,is,m,1)*f(i,it,1)
                                       -c(i,is,m,2)*f(i,it,2)
                                       -c(i,is,m,3)*f(i,it,3)
                                       -c(i,is,m,4)*f(i,it,4)
                                       -c(i,is,m,5)*f(i,it,5);
            }
        }
    }

    return;
}

} // namespace bsub_ns
