// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "vlutrp.h"
#include "vlutr.h"

namespace vlutrp_ns {

void vlutr(int& nvmax, int& n, int& nmax, int& il, int& iu, FortranArray4DRef<double> a, FortranArray4DRef<double> b, FortranArray4DRef<double> c, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    vlutr_ns::vlutr(nvmax, n, nmax, il, iu, a, b, c, nou, bou, nbuf, ibufdim);
}


void vlutrp(int& nvmax, int& n, int& nmax, int& il, int& iu, FortranArray4DRef<double> a, FortranArray4DRef<double> b, FortranArray4DRef<double> c, FortranArray4DRef<double> g, FortranArray4DRef<double> h)
{
    int k, l, izz, il1, i, ir, it, iu2, m;

    // Initialize g(izz,il,k,l) = a(izz,il,k,l) and h(izz,il,k,l) = c(izz,iu,k,l)
    for (k = 1; k <= 5; k++) {
        for (l = 1; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                g(izz,il,k,l) = a(izz,il,k,l);
                h(izz,il,k,l) = c(izz,iu,k,l);
            }
        }
    }
    il1 = il + 1;
    i   = il;

    // l-u decomposition at i=il (first part)
    for (izz = 1; izz <= n; izz++) {
        b(izz,i,1,1) = 1.e0/b(izz,i,1,1);
        b(izz,i,1,2) = b(izz,i,1,1)*b(izz,i,1,2);
        b(izz,i,1,3) = b(izz,i,1,1)*b(izz,i,1,3);
        b(izz,i,1,4) = b(izz,i,1,1)*b(izz,i,1,4);
        b(izz,i,1,5) = b(izz,i,1,1)*b(izz,i,1,5);
        b(izz,i,2,2) = 1.e0/(b(izz,i,2,2)-b(izz,i,2,1)*b(izz,i,1,2));
        b(izz,i,2,3) = b(izz,i,2,2)*(b(izz,i,2,3)-b(izz,i,2,1)*b(izz,i,1,3));
        b(izz,i,2,4) = b(izz,i,2,2)*(b(izz,i,2,4)-b(izz,i,2,1)*b(izz,i,1,4));
        b(izz,i,2,5) = b(izz,i,2,2)*(b(izz,i,2,5)-b(izz,i,2,1)*b(izz,i,1,5));
        b(izz,i,3,2) = b(izz,i,3,2)-b(izz,i,3,1)*b(izz,i,1,2);
        b(izz,i,3,3) = 1.e0/(b(izz,i,3,3)-b(izz,i,3,1)*b(izz,i,1,3)
                                         -b(izz,i,3,2)*b(izz,i,2,3));
        b(izz,i,3,4) = b(izz,i,3,3)*(b(izz,i,3,4)
                      -b(izz,i,3,1)*b(izz,i,1,4)-b(izz,i,3,2)*b(izz,i,2,4));
        b(izz,i,3,5) = b(izz,i,3,3)*(b(izz,i,3,5)
                      -b(izz,i,3,1)*b(izz,i,1,5)-b(izz,i,3,2)*b(izz,i,2,5));
    }
    for (izz = 1; izz <= n; izz++) {
        b(izz,i,4,2) = b(izz,i,4,2)-b(izz,i,4,1)*b(izz,i,1,2);
        b(izz,i,4,3) = b(izz,i,4,3)-b(izz,i,4,1)*b(izz,i,1,3)
                                   -b(izz,i,4,2)*b(izz,i,2,3);
        b(izz,i,4,4) = 1.e0/(b(izz,i,4,4)-b(izz,i,4,1)*b(izz,i,1,4)
                                         -b(izz,i,4,2)*b(izz,i,2,4)
                                         -b(izz,i,4,3)*b(izz,i,3,4));
        b(izz,i,4,5) = b(izz,i,4,4)*(b(izz,i,4,5)
                      -b(izz,i,4,1)*b(izz,i,1,5)-b(izz,i,4,2)*b(izz,i,2,5)
                      -b(izz,i,4,3)*b(izz,i,3,5));
        b(izz,i,5,2) = b(izz,i,5,2)-b(izz,i,5,1)*b(izz,i,1,2);
        b(izz,i,5,3) = b(izz,i,5,3)-b(izz,i,5,1)*b(izz,i,1,3)
                                   -b(izz,i,5,2)*b(izz,i,2,3);
        b(izz,i,5,4) = b(izz,i,5,4)-b(izz,i,5,1)*b(izz,i,1,4)
                                   -b(izz,i,5,2)*b(izz,i,2,4)
                                   -b(izz,i,5,3)*b(izz,i,3,4);
        b(izz,i,5,5) = 1.e0/(b(izz,i,5,5)-b(izz,i,5,1)*b(izz,i,1,5)
                                         -b(izz,i,5,2)*b(izz,i,2,5)
                                         -b(izz,i,5,3)*b(izz,i,3,5)
                                         -b(izz,i,5,4)*b(izz,i,4,5));
    }
    if (i == iu) goto label_1030;

    // c=ainv*c, g=ainv*g at i=il
    for (m = 1; m <= 5; m++) {
        for (izz = 1; izz <= n; izz++) {
            c(izz,i,1,m) = b(izz,i,1,1)*(c(izz,i,1,m));
            c(izz,i,2,m) = b(izz,i,2,2)*(c(izz,i,2,m)-b(izz,i,2,1)*c(izz,i,1,m));
            c(izz,i,3,m) = b(izz,i,3,3)*(c(izz,i,3,m)
                          -b(izz,i,3,1)*c(izz,i,1,m)-b(izz,i,3,2)*c(izz,i,2,m));
            c(izz,i,4,m) = b(izz,i,4,4)*(c(izz,i,4,m)
                          -b(izz,i,4,1)*c(izz,i,1,m)-b(izz,i,4,2)*c(izz,i,2,m)
                          -b(izz,i,4,3)*c(izz,i,3,m));
            c(izz,i,5,m) = b(izz,i,5,5)*(c(izz,i,5,m)
                          -b(izz,i,5,1)*c(izz,i,1,m)-b(izz,i,5,2)*c(izz,i,2,m)
                          -b(izz,i,5,3)*c(izz,i,3,m)-b(izz,i,5,4)*c(izz,i,4,m));
            c(izz,i,4,m) = c(izz,i,4,m)-b(izz,i,4,5)*c(izz,i,5,m);
            c(izz,i,3,m) = c(izz,i,3,m)-b(izz,i,3,5)*c(izz,i,5,m)
                                       -b(izz,i,3,4)*c(izz,i,4,m);
            c(izz,i,2,m) = c(izz,i,2,m)-b(izz,i,2,5)*c(izz,i,5,m)
                                       -b(izz,i,2,4)*c(izz,i,4,m)
                                       -b(izz,i,2,3)*c(izz,i,3,m);
            c(izz,i,1,m) = c(izz,i,1,m)-b(izz,i,1,5)*c(izz,i,5,m)
                                       -b(izz,i,1,4)*c(izz,i,4,m)
                                       -b(izz,i,1,3)*c(izz,i,3,m)
                                       -b(izz,i,1,2)*c(izz,i,2,m);
        }
        for (izz = 1; izz <= n; izz++) {
            g(izz,i,1,m) = b(izz,i,1,1)*(g(izz,i,1,m));
            g(izz,i,2,m) = b(izz,i,2,2)*(g(izz,i,2,m)-b(izz,i,2,1)*g(izz,i,1,m));
            g(izz,i,3,m) = b(izz,i,3,3)*(g(izz,i,3,m)
                          -b(izz,i,3,1)*g(izz,i,1,m)-b(izz,i,3,2)*g(izz,i,2,m));
            g(izz,i,4,m) = b(izz,i,4,4)*(g(izz,i,4,m)
                          -b(izz,i,4,1)*g(izz,i,1,m)-b(izz,i,4,2)*g(izz,i,2,m)
                          -b(izz,i,4,3)*g(izz,i,3,m));
            g(izz,i,5,m) = b(izz,i,5,5)*(g(izz,i,5,m)
                          -b(izz,i,5,1)*g(izz,i,1,m)-b(izz,i,5,2)*g(izz,i,2,m)
                          -b(izz,i,5,3)*g(izz,i,3,m)-b(izz,i,5,4)*g(izz,i,4,m));
            g(izz,i,4,m) = g(izz,i,4,m)-b(izz,i,4,5)*g(izz,i,5,m);
            g(izz,i,3,m) = g(izz,i,3,m)-b(izz,i,3,5)*g(izz,i,5,m)
                                       -b(izz,i,3,4)*g(izz,i,4,m);
            g(izz,i,2,m) = g(izz,i,2,m)-b(izz,i,2,5)*g(izz,i,5,m)
                                       -b(izz,i,2,4)*g(izz,i,4,m)
                                       -b(izz,i,2,3)*g(izz,i,3,m);
            g(izz,i,1,m) = g(izz,i,1,m)-b(izz,i,1,5)*g(izz,i,5,m)
                                       -b(izz,i,1,4)*g(izz,i,4,m)
                                       -b(izz,i,1,3)*g(izz,i,3,m)
                                       -b(izz,i,1,2)*g(izz,i,2,m);
        }
    }
    label_1030:;


    // forward sweep
    iu2 = iu - 2;
    if (iu2 == il) goto label_101;
    for (i = il1; i <= iu2; i++) {
        ir = i - 1;
        it = i + 1;
        // first row reduction
        for (m = 1; m <= 5; m++) {
            for (l = 1; l <= 5; l++) {
                for (izz = 1; izz <= n; izz++) {
                    b(izz,i,m,l)  = b(izz,i,m,l)
                                   -a(izz,i,m,1)*c(izz,ir,1,l)-a(izz,i,m,2)*c(izz,ir,2,l)
                                   -a(izz,i,m,3)*c(izz,ir,3,l)-a(izz,i,m,4)*c(izz,ir,4,l)
                                   -a(izz,i,m,5)*c(izz,ir,5,l);
                    g(izz,i,m,l)  = -a(izz,i,m,1)*g(izz,ir,1,l)-a(izz,i,m,2)*g(izz,ir,2,l)
                                    -a(izz,i,m,3)*g(izz,ir,3,l)-a(izz,i,m,4)*g(izz,ir,4,l)
                                    -a(izz,i,m,5)*g(izz,ir,5,l);
                    h(izz,i,m,l)  = -h(izz,ir,m,1)*c(izz,ir,1,l)-h(izz,ir,m,2)*c(izz,ir,2,l)
                                    -h(izz,ir,m,3)*c(izz,ir,3,l)-h(izz,ir,m,4)*c(izz,ir,4,l)
                                    -h(izz,ir,m,5)*c(izz,ir,5,l);
                    b(izz,iu,m,l) = b(izz,iu,m,l)
                                   -h(izz,ir,m,1)*g(izz,ir,1,l)-h(izz,ir,m,2)*g(izz,ir,2,l)
                                   -h(izz,ir,m,3)*g(izz,ir,3,l)-h(izz,ir,m,4)*g(izz,ir,4,l)
                                   -h(izz,ir,m,5)*g(izz,ir,5,l);
                }
            }
        }
        // l-u decomposition
        for (izz = 1; izz <= n; izz++) {
            b(izz,i,1,1) = 1.e0/b(izz,i,1,1);
            b(izz,i,1,2) = b(izz,i,1,1)*b(izz,i,1,2);
            b(izz,i,1,3) = b(izz,i,1,1)*b(izz,i,1,3);
            b(izz,i,1,4) = b(izz,i,1,1)*b(izz,i,1,4);
            b(izz,i,1,5) = b(izz,i,1,1)*b(izz,i,1,5);
            b(izz,i,2,2) = 1.e0/(b(izz,i,2,2)-b(izz,i,2,1)*b(izz,i,1,2));
            b(izz,i,2,3) = b(izz,i,2,2)*(b(izz,i,2,3)-b(izz,i,2,1)*b(izz,i,1,3));
            b(izz,i,2,4) = b(izz,i,2,2)*(b(izz,i,2,4)-b(izz,i,2,1)*b(izz,i,1,4));
            b(izz,i,2,5) = b(izz,i,2,2)*(b(izz,i,2,5)-b(izz,i,2,1)*b(izz,i,1,5));
            b(izz,i,3,2) = b(izz,i,3,2)-b(izz,i,3,1)*b(izz,i,1,2);
            b(izz,i,3,3) = 1.e0/(b(izz,i,3,3)-b(izz,i,3,1)*b(izz,i,1,3)
                                             -b(izz,i,3,2)*b(izz,i,2,3));
            b(izz,i,3,4) = b(izz,i,3,3)*(b(izz,i,3,4)
                          -b(izz,i,3,1)*b(izz,i,1,4)-b(izz,i,3,2)*b(izz,i,2,4));
            b(izz,i,3,5) = b(izz,i,3,3)*(b(izz,i,3,5)
                          -b(izz,i,3,1)*b(izz,i,1,5)-b(izz,i,3,2)*b(izz,i,2,5));
        }
        for (izz = 1; izz <= n; izz++) {
            b(izz,i,4,2) = b(izz,i,4,2)-b(izz,i,4,1)*b(izz,i,1,2);
            b(izz,i,4,3) = b(izz,i,4,3)-b(izz,i,4,1)*b(izz,i,1,3)
                                       -b(izz,i,4,2)*b(izz,i,2,3);
            b(izz,i,4,4) = 1.e0/(b(izz,i,4,4)-b(izz,i,4,1)*b(izz,i,1,4)
                                             -b(izz,i,4,2)*b(izz,i,2,4)
                                             -b(izz,i,4,3)*b(izz,i,3,4));
            b(izz,i,4,5) = b(izz,i,4,4)*(b(izz,i,4,5)
                          -b(izz,i,4,1)*b(izz,i,1,5)-b(izz,i,4,2)*b(izz,i,2,5)
                          -b(izz,i,4,3)*b(izz,i,3,5));
            b(izz,i,5,2) = b(izz,i,5,2)-b(izz,i,5,1)*b(izz,i,1,2);
            b(izz,i,5,3) = b(izz,i,5,3)-b(izz,i,5,1)*b(izz,i,1,3)
                                       -b(izz,i,5,2)*b(izz,i,2,3);
            b(izz,i,5,4) = b(izz,i,5,4)-b(izz,i,5,1)*b(izz,i,1,4)
                                       -b(izz,i,5,2)*b(izz,i,2,4)
                                       -b(izz,i,5,3)*b(izz,i,3,4);
            b(izz,i,5,5) = 1.e0/(b(izz,i,5,5)-b(izz,i,5,1)*b(izz,i,1,5)
                                             -b(izz,i,5,2)*b(izz,i,2,5)
                                             -b(izz,i,5,3)*b(izz,i,3,5)
                                             -b(izz,i,5,4)*b(izz,i,4,5));
        }
        // c=ainv*c
        for (m = 1; m <= 5; m++) {
            for (izz = 1; izz <= n; izz++) {
                c(izz,i,1,m) = b(izz,i,1,1)*(c(izz,i,1,m));
                c(izz,i,2,m) = b(izz,i,2,2)*(c(izz,i,2,m)-b(izz,i,2,1)*c(izz,i,1,m));
                c(izz,i,3,m) = b(izz,i,3,3)*(c(izz,i,3,m)
                              -b(izz,i,3,1)*c(izz,i,1,m)-b(izz,i,3,2)*c(izz,i,2,m));
                c(izz,i,4,m) = b(izz,i,4,4)*(c(izz,i,4,m)
                              -b(izz,i,4,1)*c(izz,i,1,m)-b(izz,i,4,2)*c(izz,i,2,m)
                              -b(izz,i,4,3)*c(izz,i,3,m));
                c(izz,i,5,m) = b(izz,i,5,5)*(c(izz,i,5,m)
                              -b(izz,i,5,1)*c(izz,i,1,m)-b(izz,i,5,2)*c(izz,i,2,m)
                              -b(izz,i,5,3)*c(izz,i,3,m)-b(izz,i,5,4)*c(izz,i,4,m));
                c(izz,i,4,m) = c(izz,i,4,m)-b(izz,i,4,5)*c(izz,i,5,m);
                c(izz,i,3,m) = c(izz,i,3,m)-b(izz,i,3,5)*c(izz,i,5,m)
                                           -b(izz,i,3,4)*c(izz,i,4,m);
                c(izz,i,2,m) = c(izz,i,2,m)-b(izz,i,2,5)*c(izz,i,5,m)
                                           -b(izz,i,2,4)*c(izz,i,4,m)
                                           -b(izz,i,2,3)*c(izz,i,3,m);
                c(izz,i,1,m) = c(izz,i,1,m)-b(izz,i,1,5)*c(izz,i,5,m)
                                           -b(izz,i,1,4)*c(izz,i,4,m)
                                           -b(izz,i,1,3)*c(izz,i,3,m)
                                           -b(izz,i,1,2)*c(izz,i,2,m);
            }
            for (izz = 1; izz <= n; izz++) {
                g(izz,i,1,m) = b(izz,i,1,1)*(g(izz,i,1,m));
                g(izz,i,2,m) = b(izz,i,2,2)*(g(izz,i,2,m)-b(izz,i,2,1)*g(izz,i,1,m));
                g(izz,i,3,m) = b(izz,i,3,3)*(g(izz,i,3,m)
                              -b(izz,i,3,1)*g(izz,i,1,m)-b(izz,i,3,2)*g(izz,i,2,m));
                g(izz,i,4,m) = b(izz,i,4,4)*(g(izz,i,4,m)
                              -b(izz,i,4,1)*g(izz,i,1,m)-b(izz,i,4,2)*g(izz,i,2,m)
                              -b(izz,i,4,3)*g(izz,i,3,m));
                g(izz,i,5,m) = b(izz,i,5,5)*(g(izz,i,5,m)
                              -b(izz,i,5,1)*g(izz,i,1,m)-b(izz,i,5,2)*g(izz,i,2,m)
                              -b(izz,i,5,3)*g(izz,i,3,m)-b(izz,i,5,4)*g(izz,i,4,m));
                g(izz,i,4,m) = g(izz,i,4,m)-b(izz,i,4,5)*g(izz,i,5,m);
                g(izz,i,3,m) = g(izz,i,3,m)-b(izz,i,3,5)*g(izz,i,5,m)
                                           -b(izz,i,3,4)*g(izz,i,4,m);
                g(izz,i,2,m) = g(izz,i,2,m)-b(izz,i,2,5)*g(izz,i,5,m)
                                           -b(izz,i,2,4)*g(izz,i,4,m)
                                           -b(izz,i,2,3)*g(izz,i,3,m);
                g(izz,i,1,m) = g(izz,i,1,m)-b(izz,i,1,5)*g(izz,i,5,m)
                                           -b(izz,i,1,4)*g(izz,i,4,m)
                                           -b(izz,i,1,3)*g(izz,i,3,m)
                                           -b(izz,i,1,2)*g(izz,i,2,m);
            }
        }
        // end of forward sweep loop body (label_1130)
    }
    label_101:;


    i  = iu - 1;
    ir = i - 1;
    it = i + 1;
    // first row reduction for i=iu-1
    for (m = 1; m <= 5; m++) {
        for (l = 1; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                b(izz,i,m,l)  = b(izz,i,m,l)
                               -a(izz,i,m,1)*c(izz,ir,1,l)-a(izz,i,m,2)*c(izz,ir,2,l)
                               -a(izz,i,m,3)*c(izz,ir,3,l)-a(izz,i,m,4)*c(izz,ir,4,l)
                               -a(izz,i,m,5)*c(izz,ir,5,l);
                c(izz,i,m,l)  = c(izz,i,m,l)
                               -a(izz,i,m,1)*g(izz,ir,1,l)-a(izz,i,m,2)*g(izz,ir,2,l)
                               -a(izz,i,m,3)*g(izz,ir,3,l)-a(izz,i,m,4)*g(izz,ir,4,l)
                               -a(izz,i,m,5)*g(izz,ir,5,l);
                a(izz,iu,m,l) = a(izz,iu,m,l)
                               -h(izz,ir,m,1)*c(izz,ir,1,l)-h(izz,ir,m,2)*c(izz,ir,2,l)
                               -h(izz,ir,m,3)*c(izz,ir,3,l)-h(izz,ir,m,4)*c(izz,ir,4,l)
                               -h(izz,ir,m,5)*c(izz,ir,5,l);
                b(izz,iu,m,l) = b(izz,iu,m,l)
                               -h(izz,ir,m,1)*g(izz,ir,1,l)-h(izz,ir,m,2)*g(izz,ir,2,l)
                               -h(izz,ir,m,3)*g(izz,ir,3,l)-h(izz,ir,m,4)*g(izz,ir,4,l)
                               -h(izz,ir,m,5)*g(izz,ir,5,l);
            }
        }
    }
    // l-u decomposition for i=iu-1
    for (izz = 1; izz <= n; izz++) {
        b(izz,i,1,1) = 1.e0/b(izz,i,1,1);
        b(izz,i,1,2) = b(izz,i,1,1)*b(izz,i,1,2);
        b(izz,i,1,3) = b(izz,i,1,1)*b(izz,i,1,3);
        b(izz,i,1,4) = b(izz,i,1,1)*b(izz,i,1,4);
        b(izz,i,1,5) = b(izz,i,1,1)*b(izz,i,1,5);
        b(izz,i,2,2) = 1.e0/(b(izz,i,2,2)-b(izz,i,2,1)*b(izz,i,1,2));
        b(izz,i,2,3) = b(izz,i,2,2)*(b(izz,i,2,3)-b(izz,i,2,1)*b(izz,i,1,3));
        b(izz,i,2,4) = b(izz,i,2,2)*(b(izz,i,2,4)-b(izz,i,2,1)*b(izz,i,1,4));
        b(izz,i,2,5) = b(izz,i,2,2)*(b(izz,i,2,5)-b(izz,i,2,1)*b(izz,i,1,5));
        b(izz,i,3,2) = b(izz,i,3,2)-b(izz,i,3,1)*b(izz,i,1,2);
        b(izz,i,3,3) = 1.e0/(b(izz,i,3,3)-b(izz,i,3,1)*b(izz,i,1,3)
                                         -b(izz,i,3,2)*b(izz,i,2,3));
        b(izz,i,3,4) = b(izz,i,3,3)*(b(izz,i,3,4)
                      -b(izz,i,3,1)*b(izz,i,1,4)-b(izz,i,3,2)*b(izz,i,2,4));
        b(izz,i,3,5) = b(izz,i,3,3)*(b(izz,i,3,5)
                      -b(izz,i,3,1)*b(izz,i,1,5)-b(izz,i,3,2)*b(izz,i,2,5));
    }
    for (izz = 1; izz <= n; izz++) {
        b(izz,i,4,2) = b(izz,i,4,2)-b(izz,i,4,1)*b(izz,i,1,2);
        b(izz,i,4,3) = b(izz,i,4,3)-b(izz,i,4,1)*b(izz,i,1,3)
                                   -b(izz,i,4,2)*b(izz,i,2,3);
        b(izz,i,4,4) = 1.e0/(b(izz,i,4,4)-b(izz,i,4,1)*b(izz,i,1,4)
                                         -b(izz,i,4,2)*b(izz,i,2,4)
                                         -b(izz,i,4,3)*b(izz,i,3,4));
        b(izz,i,4,5) = b(izz,i,4,4)*(b(izz,i,4,5)
                      -b(izz,i,4,1)*b(izz,i,1,5)-b(izz,i,4,2)*b(izz,i,2,5)
                      -b(izz,i,4,3)*b(izz,i,3,5));
        b(izz,i,5,2) = b(izz,i,5,2)-b(izz,i,5,1)*b(izz,i,1,2);
        b(izz,i,5,3) = b(izz,i,5,3)-b(izz,i,5,1)*b(izz,i,1,3)
                                   -b(izz,i,5,2)*b(izz,i,2,3);
        b(izz,i,5,4) = b(izz,i,5,4)-b(izz,i,5,1)*b(izz,i,1,4)
                                   -b(izz,i,5,2)*b(izz,i,2,4)
                                   -b(izz,i,5,3)*b(izz,i,3,4);
        b(izz,i,5,5) = 1.e0/(b(izz,i,5,5)-b(izz,i,5,1)*b(izz,i,1,5)
                                         -b(izz,i,5,2)*b(izz,i,2,5)
                                         -b(izz,i,5,3)*b(izz,i,3,5)
                                         -b(izz,i,5,4)*b(izz,i,4,5));
    }
    // c=ainv*c for i=iu-1
    for (m = 1; m <= 5; m++) {
        for (izz = 1; izz <= n; izz++) {
            c(izz,i,1,m) = b(izz,i,1,1)*(c(izz,i,1,m));
            c(izz,i,2,m) = b(izz,i,2,2)*(c(izz,i,2,m)-b(izz,i,2,1)*c(izz,i,1,m));
            c(izz,i,3,m) = b(izz,i,3,3)*(c(izz,i,3,m)
                          -b(izz,i,3,1)*c(izz,i,1,m)-b(izz,i,3,2)*c(izz,i,2,m));
            c(izz,i,4,m) = b(izz,i,4,4)*(c(izz,i,4,m)
                          -b(izz,i,4,1)*c(izz,i,1,m)-b(izz,i,4,2)*c(izz,i,2,m)
                          -b(izz,i,4,3)*c(izz,i,3,m));
            c(izz,i,5,m) = b(izz,i,5,5)*(c(izz,i,5,m)
                          -b(izz,i,5,1)*c(izz,i,1,m)-b(izz,i,5,2)*c(izz,i,2,m)
                          -b(izz,i,5,3)*c(izz,i,3,m)-b(izz,i,5,4)*c(izz,i,4,m));
            c(izz,i,4,m) = c(izz,i,4,m)-b(izz,i,4,5)*c(izz,i,5,m);
            c(izz,i,3,m) = c(izz,i,3,m)-b(izz,i,3,5)*c(izz,i,5,m)
                                       -b(izz,i,3,4)*c(izz,i,4,m);
            c(izz,i,2,m) = c(izz,i,2,m)-b(izz,i,2,5)*c(izz,i,5,m)
                                       -b(izz,i,2,4)*c(izz,i,4,m)
                                       -b(izz,i,2,3)*c(izz,i,3,m);
            c(izz,i,1,m) = c(izz,i,1,m)-b(izz,i,1,5)*c(izz,i,5,m)
                                       -b(izz,i,1,4)*c(izz,i,4,m)
                                       -b(izz,i,1,3)*c(izz,i,3,m)
                                       -b(izz,i,1,2)*c(izz,i,2,m);
        }
    }


    i  = iu;
    ir = i - 1;
    it = i + 1;
    // row reduction for i=iu
    for (m = 1; m <= 5; m++) {
        for (l = 1; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                b(izz,i,m,l) = b(izz,i,m,l)
                              -a(izz,i,m,1)*c(izz,ir,1,l)-a(izz,i,m,2)*c(izz,ir,2,l)
                              -a(izz,i,m,3)*c(izz,ir,3,l)-a(izz,i,m,4)*c(izz,ir,4,l)
                              -a(izz,i,m,5)*c(izz,ir,5,l);
            }
        }
    }
    // l-u decomposition for i=iu
    for (izz = 1; izz <= n; izz++) {
        b(izz,i,1,1) = 1.e0/b(izz,i,1,1);
        b(izz,i,1,2) = b(izz,i,1,1)*b(izz,i,1,2);
        b(izz,i,1,3) = b(izz,i,1,1)*b(izz,i,1,3);
        b(izz,i,1,4) = b(izz,i,1,1)*b(izz,i,1,4);
        b(izz,i,1,5) = b(izz,i,1,1)*b(izz,i,1,5);
        b(izz,i,2,2) = 1.e0/(b(izz,i,2,2)-b(izz,i,2,1)*b(izz,i,1,2));
        b(izz,i,2,3) = b(izz,i,2,2)*(b(izz,i,2,3)-b(izz,i,2,1)*b(izz,i,1,3));
        b(izz,i,2,4) = b(izz,i,2,2)*(b(izz,i,2,4)-b(izz,i,2,1)*b(izz,i,1,4));
        b(izz,i,2,5) = b(izz,i,2,2)*(b(izz,i,2,5)-b(izz,i,2,1)*b(izz,i,1,5));
        b(izz,i,3,2) = b(izz,i,3,2)-b(izz,i,3,1)*b(izz,i,1,2);
        b(izz,i,3,3) = 1.e0/(b(izz,i,3,3)-b(izz,i,3,1)*b(izz,i,1,3)
                                         -b(izz,i,3,2)*b(izz,i,2,3));
        b(izz,i,3,4) = b(izz,i,3,3)*(b(izz,i,3,4)
                      -b(izz,i,3,1)*b(izz,i,1,4)-b(izz,i,3,2)*b(izz,i,2,4));
        b(izz,i,3,5) = b(izz,i,3,3)*(b(izz,i,3,5)
                      -b(izz,i,3,1)*b(izz,i,1,5)-b(izz,i,3,2)*b(izz,i,2,5));
    }
    for (izz = 1; izz <= n; izz++) {
        b(izz,i,4,2) = b(izz,i,4,2)-b(izz,i,4,1)*b(izz,i,1,2);
        b(izz,i,4,3) = b(izz,i,4,3)-b(izz,i,4,1)*b(izz,i,1,3)
                                   -b(izz,i,4,2)*b(izz,i,2,3);
        b(izz,i,4,4) = 1.e0/(b(izz,i,4,4)-b(izz,i,4,1)*b(izz,i,1,4)
                                         -b(izz,i,4,2)*b(izz,i,2,4)
                                         -b(izz,i,4,3)*b(izz,i,3,4));
        b(izz,i,4,5) = b(izz,i,4,4)*(b(izz,i,4,5)
                      -b(izz,i,4,1)*b(izz,i,1,5)-b(izz,i,4,2)*b(izz,i,2,5)
                      -b(izz,i,4,3)*b(izz,i,3,5));
        b(izz,i,5,2) = b(izz,i,5,2)-b(izz,i,5,1)*b(izz,i,1,2);
        b(izz,i,5,3) = b(izz,i,5,3)-b(izz,i,5,1)*b(izz,i,1,3)
                                   -b(izz,i,5,2)*b(izz,i,2,3);
        b(izz,i,5,4) = b(izz,i,5,4)-b(izz,i,5,1)*b(izz,i,1,4)
                                   -b(izz,i,5,2)*b(izz,i,2,4)
                                   -b(izz,i,5,3)*b(izz,i,3,4);
        b(izz,i,5,5) = 1.e0/(b(izz,i,5,5)-b(izz,i,5,1)*b(izz,i,1,5)
                                         -b(izz,i,5,2)*b(izz,i,2,5)
                                         -b(izz,i,5,3)*b(izz,i,3,5)
                                         -b(izz,i,5,4)*b(izz,i,4,5));
    }
    return;
}

} // namespace vlutrp_ns
