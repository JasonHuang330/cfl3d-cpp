// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "tmetric.h"
#include <cstdio>

namespace tmetric_ns {

void tmetric(int& jdim, int& kdim, int& idim,
             FortranArray3DRef<double> sj,
             FortranArray3DRef<double> sk,
             FortranArray3DRef<double> si,
             FortranArray3DRef<double> x,
             FortranArray3DRef<double> y,
             FortranArray3DRef<double> z,
             FortranArray2DRef<double> t,
             FortranArray2DRef<double> t1,
             FortranArray2DRef<double> t2,
             FortranArray3DRef<double> t3,
             int& nbl)
{
    int kdim1, idim1;
    int i, izz, k;
    int js, jk1, jki1, jk1i1;
    int j1k, j1ki1;
    int j1k1;
    int n, jk;

    kdim1 = kdim - 1;
    idim1 = idim - 1;

    //
    // Loop over i = 1..idim1
    //
    for (i = 1; i <= idim1; i++) {

        // **********************************************
        // unsteady sj terms
        // **********************************************

        // face-average grid speeds
        js    = jdim * kdim * (i - 1) + 1;
        jk1   = js + jdim;
        jki1  = jdim * kdim * i + 1;
        jk1i1 = jki1 + jdim;
        n     = jdim * kdim1;

        for (izz = 1; izz <= n; izz++) {
            t1(izz, 1) = 0.25 * (t(izz + js - 1, 1)   + t(izz + jk1 - 1, 1)
                                + t(izz + jki1 - 1, 1) + t(izz + jk1i1 - 1, 1));
            t1(izz, 2) = 0.25 * (t(izz + js - 1, 2)   + t(izz + jk1 - 1, 2)
                                + t(izz + jki1 - 1, 2) + t(izz + jk1i1 - 1, 2));
            t1(izz, 3) = 0.25 * (t(izz + js - 1, 3)   + t(izz + jk1 - 1, 3)
                                + t(izz + jki1 - 1, 3) + t(izz + jk1i1 - 1, 3));

            sj(izz, i, 5) = - t1(izz, 1) * sj(izz, i, 1)
                            - t1(izz, 2) * sj(izz, i, 2)
                            - t1(izz, 3) * sj(izz, i, 3);
        }

        // zero out sj(j,k=kdim,i,5)
        for (izz = 1; izz <= jdim; izz++) {
            sj(izz + n, i, 5) = 0.0e0;
        }

        // **********************************************
        // unsteady sk terms
        // **********************************************

        // face-average grid speeds
        j1k   = js + 1;
        j1ki1 = jki1 + 1;
        n     = jdim * kdim - 1;

        for (izz = 1; izz <= n; izz++) {
            t1(izz, 1) = 0.25 * (t(izz + js - 1, 1)   + t(izz + j1k - 1, 1)
                                + t(izz + jki1 - 1, 1) + t(izz + j1ki1 - 1, 1));
            t1(izz, 2) = 0.25 * (t(izz + js - 1, 2)   + t(izz + j1k - 1, 2)
                                + t(izz + jki1 - 1, 2) + t(izz + j1ki1 - 1, 2));
            t1(izz, 3) = 0.25 * (t(izz + js - 1, 3)   + t(izz + j1k - 1, 3)
                                + t(izz + jki1 - 1, 3) + t(izz + j1ki1 - 1, 3));
        }

        n = jdim * kdim;
        for (izz = 1; izz <= n; izz++) {
            sk(izz, i, 5) = - t1(izz, 1) * sk(izz, i, 1)
                            - t1(izz, 2) * sk(izz, i, 2)
                            - t1(izz, 3) * sk(izz, i, 3);
        }

        // zero out sk(j=jdim,k,i,5)
        for (k = 1; k <= kdim; k++) {
            jk = jdim * k;
            sk(jk, i, 5) = 0.0;
        }

    } // end do 50

    //
    // Loop over i = 1..idim
    //
    for (i = 1; i <= idim; i++) {

        // **********************************************
        // unsteady si terms
        // **********************************************

        // face-average grid speeds
        js   = jdim * kdim * (i - 1) + 1;
        j1k  = js + 1;
        jk1  = js + jdim;
        j1k1 = jk1 + 1;
        n    = jdim * (kdim - 1) - 1;

        for (izz = 1; izz <= n; izz++) {
            t1(izz, 1) = 0.25 * (t(izz + js - 1, 1)  + t(izz + j1k - 1, 1)
                                + t(izz + jk1 - 1, 1) + t(izz + j1k1 - 1, 1));
            t1(izz, 2) = 0.25 * (t(izz + js - 1, 2)  + t(izz + j1k - 1, 2)
                                + t(izz + jk1 - 1, 2) + t(izz + j1k1 - 1, 2));
            t1(izz, 3) = 0.25 * (t(izz + js - 1, 3)  + t(izz + j1k - 1, 3)
                                + t(izz + jk1 - 1, 3) + t(izz + j1k1 - 1, 3));
        }

        n = jdim * kdim1;
        for (izz = 1; izz <= n; izz++) {
            si(izz, i, 5) = - t1(izz, 1) * si(izz, i, 1)
                            - t1(izz, 2) * si(izz, i, 2)
                            - t1(izz, 3) * si(izz, i, 3);
        }

        // zero out si(j,k=kdim,i,5)
        for (izz = 1; izz <= jdim; izz++) {
            si(izz + n, i, 5) = 0.0e0;
        }

        // zero out si(j=jdim,k,i,5)
        for (k = 1; k <= kdim1; k++) {
            jk = jdim * k;
            si(jk, i, 5) = 0.0e0;
        }

    } // end do 40

    return;
}

} // namespace tmetric_ns
