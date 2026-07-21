// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "vlutr.h"

namespace vlutr_ns {

void vlutr(int& nvmax, int& n, int& nmax, int& il, int& iu,
           FortranArray4DRef<double> a,
           FortranArray4DRef<double> b,
           FortranArray4DRef<double> c,
           FortranArray1DRef<int> nou,
           FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim)
{
    // LU decomposition for a block 5x5 tridiagonal system of equations.
    // a, b, c are (nvmax, nmax, 5, 5) arrays
    // il and iu are starting and finishing indices
    // solution is by upper triangularization with unity diagonal
    // block inversions use nonpivoted LU decomposition

    int il1 = il + 1;
    int i   = il;

    // L-U decomposition of first block
    // loop 1000
    for (int izz = 1; izz <= n; izz++) {
        b(izz,i,1,1) = 1.0/b(izz,i,1,1);
        b(izz,i,1,2) = b(izz,i,1,1)*b(izz,i,1,2);
        b(izz,i,1,3) = b(izz,i,1,1)*b(izz,i,1,3);
        b(izz,i,1,4) = b(izz,i,1,1)*b(izz,i,1,4);
        b(izz,i,1,5) = b(izz,i,1,1)*b(izz,i,1,5);
        // b(izz,i,2,1) = b(izz,i,2,1)
        b(izz,i,2,2) = 1.0/(b(izz,i,2,2)-b(izz,i,2,1)*b(izz,i,1,2));
        b(izz,i,2,3) = b(izz,i,2,2)*(b(izz,i,2,3)
                      -b(izz,i,2,1)* b(izz,i,1,3));
        b(izz,i,2,4) = b(izz,i,2,2)*(b(izz,i,2,4)
                      -b(izz,i,2,1)* b(izz,i,1,4));
        b(izz,i,2,5) = b(izz,i,2,2)*(b(izz,i,2,5)
                      -b(izz,i,2,1)* b(izz,i,1,5));
        // b(izz,i,3,1) = b(izz,i,3,1)
        b(izz,i,3,2) = b(izz,i,3,2)-b(izz,i,3,1)*b(izz,i,1,2);
        b(izz,i,3,3) = 1.0/(b(izz,i,3,3)-b(izz,i,3,1)*b(izz,i,1,3)
                                         -b(izz,i,3,2)*b(izz,i,2,3));
        b(izz,i,3,4) = b(izz,i,3,3)*(b(izz,i,3,4)
                      -b(izz,i,3,1)* b(izz,i,1,4)
                      -b(izz,i,3,2)* b(izz,i,2,4));
        b(izz,i,3,5) = b(izz,i,3,3)*(b(izz,i,3,5)
                      -b(izz,i,3,1)*b(izz,i,1,5)
                      -b(izz,i,3,2)*b(izz,i,2,5));
        // b(izz,i,4,1) = b(izz,i,4,1)
    } // end loop 1000

    // loop 3000
    for (int izz = 1; izz <= n; izz++) {
        b(izz,i,4,2) = b(izz,i,4,2)-b(izz,i,4,1)*b(izz,i,1,2);
        b(izz,i,4,3) = b(izz,i,4,3)-b(izz,i,4,1)*b(izz,i,1,3)
                                    -b(izz,i,4,2)*b(izz,i,2,3);
        b(izz,i,4,4) = 1.0/(b(izz,i,4,4)-b(izz,i,4,1)*b(izz,i,1,4)
                            -b(izz,i,4,2)*b(izz,i,2,4)
                            -b(izz,i,4,3)*b(izz,i,3,4));
        b(izz,i,4,5) = b(izz,i,4,4)*(b(izz,i,4,5)
                                    -b(izz,i,4,1)*b(izz,i,1,5)
                                    -b(izz,i,4,2)* b(izz,i,2,5)
                                    -b(izz,i,4,3)*b(izz,i,3,5));
        // b(izz,i,5,1) = b(izz,i,5,1)
        b(izz,i,5,2) = b(izz,i,5,2)-b(izz,i,5,1)*b(izz,i,1,2);
        b(izz,i,5,3) = b(izz,i,5,3)-b(izz,i,5,1)*b(izz,i,1,3)
                                    -b(izz,i,5,2)*b(izz,i,2,3);
        b(izz,i,5,4) = b(izz,i,5,4)-b(izz,i,5,1)*b(izz,i,1,4)
                                    -b(izz,i,5,2)*b(izz,i,2,4)
                                    -b(izz,i,5,3)*b(izz,i,3,4);
        b(izz,i,5,5) = 1.0/(b(izz,i,5,5)-b(izz,i,5,1)*b(izz,i,1,5)
                                         -b(izz,i,5,2)*b(izz,i,2,5)
                                         -b(izz,i,5,3)*b(izz,i,3,5)
                                         -b(izz,i,5,4)*b(izz,i,4,5));
    } // end loop 3000

    if (i != iu) {
        // c = ainv*c  (loops 1001, 1010)
        for (int m = 1; m <= 5; m++) {
            for (int izz = 1; izz <= n; izz++) {
                c(izz,i,1,m) = b(izz,i,1,1)*(c(izz,i,1,m));
                c(izz,i,2,m) = b(izz,i,2,2)*(c(izz,i,2,m)
                              -b(izz,i,2,1)* c(izz,i,1,m));
                c(izz,i,3,m) = b(izz,i,3,3)*(c(izz,i,3,m)
                              -b(izz,i,3,1)* c(izz,i,1,m)
                              -b(izz,i,3,2)* c(izz,i,2,m));
                c(izz,i,4,m) = b(izz,i,4,4)*(c(izz,i,4,m)
                              -b(izz,i,4,1)* c(izz,i,1,m)
                              -b(izz,i,4,2)* c(izz,i,2,m)
                              -b(izz,i,4,3)* c(izz,i,3,m));
                c(izz,i,5,m) = b(izz,i,5,5)*(c(izz,i,5,m)
                              -b(izz,i,5,1)* c(izz,i,1,m)
                              -b(izz,i,5,2)* c(izz,i,2,m)
                              -b(izz,i,5,3)* c(izz,i,3,m)
                              -b(izz,i,5,4)* c(izz,i,4,m));
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
            } // end loop 1001
        } // end loop 1010
    } // end if (i != iu) [label 1030]

    // Forward sweep: loop 100
    for (i = il1; i <= iu; i++) {
        int ir = i - 1;
        int it = i + 1;

        // First row reduction: b(i) -= a(i)*c(ir)  (loops 20, 1002)
        for (int m = 1; m <= 5; m++) {
            for (int l = 1; l <= 5; l++) {
                for (int izz = 1; izz <= n; izz++) {
                    b(izz,i,m,l) = b(izz,i,m,l)-a(izz,i,m,1)*c(izz,ir,1,l)
                                               -a(izz,i,m,2)*c(izz,ir,2,l)
                                               -a(izz,i,m,3)*c(izz,ir,3,l)
                                               -a(izz,i,m,4)*c(izz,ir,4,l)
                                               -a(izz,i,m,5)*c(izz,ir,5,l);
                } // end loop 1002
            }
        } // end loop 20

        // L-U decomposition of current block (loop 1003)
        for (int izz = 1; izz <= n; izz++) {
            b(izz,i,1,1) = 1.0/b(izz,i,1,1);
            b(izz,i,1,2) = b(izz,i,1,1)*b(izz,i,1,2);
            b(izz,i,1,3) = b(izz,i,1,1)*b(izz,i,1,3);
            b(izz,i,1,4) = b(izz,i,1,1)*b(izz,i,1,4);
            b(izz,i,1,5) = b(izz,i,1,1)*b(izz,i,1,5);
            // b(izz,i,2,1) = b(izz,i,2,1)
            b(izz,i,2,2) = 1.0/(b(izz,i,2,2)-b(izz,i,2,1)*b(izz,i,1,2));
            b(izz,i,2,3) = b(izz,i,2,2)*(b(izz,i,2,3)
                          -b(izz,i,2,1)* b(izz,i,1,3));
            b(izz,i,2,4) = b(izz,i,2,2)*(b(izz,i,2,4)
                          -b(izz,i,2,1)* b(izz,i,1,4));
            b(izz,i,2,5) = b(izz,i,2,2)*(b(izz,i,2,5)
                          -b(izz,i,2,1)* b(izz,i,1,5));
            // b(izz,i,3,1) = b(izz,i,3,1)
            b(izz,i,3,2) = b(izz,i,3,2)-b(izz,i,3,1)*b(izz,i,1,2);
            b(izz,i,3,3) = 1.0/(b(izz,i,3,3)-b(izz,i,3,1)*b(izz,i,1,3)
                                             -b(izz,i,3,2)*b(izz,i,2,3));
            b(izz,i,3,4) = b(izz,i,3,3)*(b(izz,i,3,4)
                          -b(izz,i,3,1)* b(izz,i,1,4)
                          -b(izz,i,3,2)* b(izz,i,2,4));
            b(izz,i,3,5) = b(izz,i,3,3)*(b(izz,i,3,5)
                          -b(izz,i,3,1)* b(izz,i,1,5)
                          -b(izz,i,3,2)* b(izz,i,2,5));
            // b(izz,i,4,1) = b(izz,i,4,1)
        } // end loop 1003

        // loop 3003
        for (int izz = 1; izz <= n; izz++) {
            b(izz,i,4,2) = b(izz,i,4,2)-b(izz,i,4,1)*b(izz,i,1,2);
            b(izz,i,4,3) = b(izz,i,4,3)-b(izz,i,4,1)*b(izz,i,1,3)
                                        -b(izz,i,4,2)*b(izz,i,2,3);
            b(izz,i,4,4) = 1.0/(b(izz,i,4,4)-b(izz,i,4,1)*b(izz,i,1,4)
                                             -b(izz,i,4,2)*b(izz,i,2,4)
                                             -b(izz,i,4,3)*b(izz,i,3,4));
            b(izz,i,4,5) = b(izz,i,4,4)*(b(izz,i,4,5)
                          -b(izz,i,4,1)* b(izz,i,1,5)
                          -b(izz,i,4,2)* b(izz,i,2,5)
                          -b(izz,i,4,3)* b(izz,i,3,5));
            // b(izz,i,5,1) = b(izz,i,5,1)
            b(izz,i,5,2) = b(izz,i,5,2)-b(izz,i,5,1)*b(izz,i,1,2);
            b(izz,i,5,3) = b(izz,i,5,3)-b(izz,i,5,1)*b(izz,i,1,3)
                                        -b(izz,i,5,2)*b(izz,i,2,3);
            b(izz,i,5,4) = b(izz,i,5,4)-b(izz,i,5,1)*b(izz,i,1,4)
                                        -b(izz,i,5,2)*b(izz,i,2,4)
                                        -b(izz,i,5,3)*b(izz,i,3,4);
            b(izz,i,5,5) = 1.0/(b(izz,i,5,5)-b(izz,i,5,1)*b(izz,i,1,5)
                                             -b(izz,i,5,2)*b(izz,i,2,5)
                                             -b(izz,i,5,3)*b(izz,i,3,5)
                                             -b(izz,i,5,4)*b(izz,i,4,5));
        } // end loop 3003

        if (i != iu) {
            // c = ainv*c  (loops 1004, 1110)
            for (int m = 1; m <= 5; m++) {
                for (int izz = 1; izz <= n; izz++) {
                    c(izz,i,1,m) = b(izz,i,1,1)*(c(izz,i,1,m));
                    c(izz,i,2,m) = b(izz,i,2,2)*(c(izz,i,2,m)
                                  -b(izz,i,2,1)* c(izz,i,1,m));
                    c(izz,i,3,m) = b(izz,i,3,3)*(c(izz,i,3,m)
                                  -b(izz,i,3,1)* c(izz,i,1,m)
                                  -b(izz,i,3,2)* c(izz,i,2,m));
                    c(izz,i,4,m) = b(izz,i,4,4)*(c(izz,i,4,m)
                                  -b(izz,i,4,1)* c(izz,i,1,m)
                                  -b(izz,i,4,2)* c(izz,i,2,m)
                                  -b(izz,i,4,3)* c(izz,i,3,m));
                    c(izz,i,5,m) = b(izz,i,5,5)*(c(izz,i,5,m)
                                  -b(izz,i,5,1)* c(izz,i,1,m)
                                  -b(izz,i,5,2)* c(izz,i,2,m)
                                  -b(izz,i,5,3)* c(izz,i,3,m)
                                  -b(izz,i,5,4)*c(izz,i,4,m));
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
                } // end loop 1004
            } // end loop 1110
        } // end if (i != iu) [label 1130]
    } // end loop 100

    return;
}

} // namespace vlutr_ns
