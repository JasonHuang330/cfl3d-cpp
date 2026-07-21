// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bsubp.h"
#include "bsub.h"

namespace bsubp_ns {

void bsub(int& id1, int& id2, FortranArray4DRef<double> a, FortranArray4DRef<double> b, FortranArray4DRef<double> c, FortranArray3DRef<double> f, int& i1, int& i2, int& il, int& iu)
{
    bsub_ns::bsub(id1, id2, a, b, c, f, i1, i2, il, iu);
}

void bsubp(int& id1, int& id2, FortranArray4DRef<double> a, FortranArray4DRef<double> b, FortranArray4DRef<double> c, FortranArray3DRef<double> f, int& i1, int& i2, int& il, int& iu, FortranArray4DRef<double> g, FortranArray4DRef<double> h)
{
    int il1, is, iux, ir, it, m, i, ix, iqq, iu2, il11;
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
    iux = iu - 1;
    for (is = il1; is <= iux; is++) {
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

    is = iu;
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
    iu2 = iu - 2;
    for (ix = il; ix <= iu2; ix++) {
        for (m = 1; m <= 5; m++) {
            for (i = i1; i <= i2; i++) {
                f(i,is,m) = f(i,is,m)-h(i,ix,m,1)*f(i,ix,1)
                                      -h(i,ix,m,2)*f(i,ix,2)
                                      -h(i,ix,m,3)*f(i,ix,3)
                                      -h(i,ix,m,4)*f(i,ix,4)
                                      -h(i,ix,m,5)*f(i,ix,5);
            }
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

    // back substitution
    iux = il1;
    for (iqq = il1; iqq <= iux; iqq++) {
        is = il + iu - iqq;
        it = is + 1;
        for (m = 1; m <= 5; m++) {
            for (i = i1; i <= i2; i++) {
                f(i,is,m) = f(i,is,m)-c(i,is,m,1)*f(i,it,1)
                                      -c(i,is,m,2)*f(i,it,2)
                                      -c(i,is,m,3)*f(i,it,3)
                                      -c(i,is,m,4)*f(i,it,4)
                                      -c(i,is,m,5)*f(i,it,5);
            }
        }
    }

    il11 = il1 + 1;
    for (iqq = il11; iqq <= iu; iqq++) {
        is = il + iu - iqq;
        it = is + 1;
        for (m = 1; m <= 5; m++) {
            for (i = i1; i <= i2; i++) {
                f(i,is,m) = f(i,is,m)-c(i,is,m,1)*f(i,it,1)
                                      -c(i,is,m,2)*f(i,it,2)
                                      -c(i,is,m,3)*f(i,it,3)
                                      -c(i,is,m,4)*f(i,it,4)
                                      -c(i,is,m,5)*f(i,it,5)-g(i,is,m,1)*f(i,iu,1)
                                                             -g(i,is,m,2)*f(i,iu,2)
                                                             -g(i,is,m,3)*f(i,iu,3)
                                                             -g(i,is,m,4)*f(i,iu,4)
                                                             -g(i,is,m,5)*f(i,iu,5);
            }
        }
    }
}

} // namespace bsubp_ns
