// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "cellvol.h"
#include "termn8.h"
#include "ccomplex.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace cellvol_ns {

void cellvol(int& jdim, int& kdim, int& idim,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
             FortranArray3DRef<double> sj, FortranArray3DRef<double> sk, FortranArray3DRef<double> si,
             FortranArray3DRef<double> vol,
             FortranArray2DRef<double> t,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& myid,
             FortranArray1DRef<int> mblk2nd, int& maxbl, int& nbl,
             int& iflagv1, int& iflagv,
             int& imin, int& imax, int& jmin, int& jmax, int& kmin, int& kmax)
{
    // COMMON block references
    int32_t& negvol = cmn_deformz.negvol;
    int32_t& iexp   = cmn_zero.iexp;

    double eps = std::pow(10.0, (double)(-iexp + 1));
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    //
    // directed area and magnitude (i=constant face)
    //
    int i = 1;
    int n = jdim * kdim1 - 1;

    for (int izz = 1; izz <= n; izz++) {
        t(izz, 7) = si(izz, i, 1) * si(izz, i, 4);
        t(izz, 8) = si(izz, i, 2) * si(izz, i, 4);
        t(izz, 9) = si(izz, i, 3) * si(izz, i, 4);
        //
        // average point in i=constant face
        //
        t(izz, 13) = x(izz, 1, i) + x(izz+1, 1, i) + x(izz, 2, i) + x(izz+1, 2, i);
        t(izz, 14) = y(izz, 1, i) + y(izz+1, 1, i) + y(izz, 2, i) + y(izz+1, 2, i);
        t(izz, 15) = z(izz, 1, i) + z(izz+1, 1, i) + z(izz, 2, i) + z(izz+1, 2, i);
    }

    double t1 = -1.e0 / 24.e0;

    for (i = 1; i <= idim1; i++) {
        //
        // average point in volume
        //
        int l = i + 1;
        n = jdim * kdim1 - 1;

        for (int izz = 1; izz <= n; izz++) {
            t(izz, 10) = x(izz, 1, l) + x(izz+1, 1, l) + x(izz, 2, l) + x(izz+1, 2, l)
                       + t(izz, 13);
            t(izz, 11) = y(izz, 1, l) + y(izz+1, 1, l) + y(izz, 2, l) + y(izz+1, 2, l)
                       + t(izz, 14);
            t(izz, 12) = z(izz, 1, l) + z(izz+1, 1, l) + z(izz, 2, l) + z(izz+1, 2, l)
                       + t(izz, 15);
            //
            // accumulate volume as sum of pentahedrons
            //
            vol(izz, 1, i) = (2.e0 * t(izz, 13) - t(izz, 10)) * t(izz, 7)
                           + (2.e0 * t(izz, 14) - t(izz, 11)) * t(izz, 8)
                           + (2.e0 * t(izz, 15) - t(izz, 12)) * t(izz, 9);
        }

        //
        // directed area (j=constant face)
        //
        n = jdim * kdim1;

        for (int izz = 1; izz <= n; izz++) {
            t(izz, 7) = sj(izz, i, 1) * sj(izz, i, 4);
            t(izz, 8) = sj(izz, i, 2) * sj(izz, i, 4);
            t(izz, 9) = sj(izz, i, 3) * sj(izz, i, 4);
            //
            // average point in j=constant face
            //
            t(izz, 4) = x(izz, 1, i) + x(izz, 1, i+1) + x(izz, 2, i) + x(izz, 2, i+1);
            t(izz, 5) = y(izz, 1, i) + y(izz, 1, i+1) + y(izz, 2, i) + y(izz, 2, i+1);
            t(izz, 6) = z(izz, 1, i) + z(izz, 1, i+1) + z(izz, 2, i) + z(izz, 2, i+1);
        }

        n = n - 1;

        for (int izz = 1; izz <= n; izz++) {
            vol(izz, 1, i) = vol(izz, 1, i)
                           + (2.e0 * t(izz,   4) - t(izz, 10)) * t(izz,   7)
                           + (2.e0 * t(izz,   5) - t(izz, 11)) * t(izz,   8)
                           + (2.e0 * t(izz,   6) - t(izz, 12)) * t(izz,   9)
                           - (2.e0 * t(izz+1, 4) - t(izz, 10)) * t(izz+1, 7)
                           - (2.e0 * t(izz+1, 5) - t(izz, 11)) * t(izz+1, 8)
                           - (2.e0 * t(izz+1, 6) - t(izz, 12)) * t(izz+1, 9);
        }

        //
        // directed area and magnitude (k=constant face)
        //
        n = jdim * kdim - 1;

        for (int izz = 1; izz <= n; izz++) {
            t(izz, 7) = sk(izz, i, 1) * sk(izz, i, 4);
            t(izz, 8) = sk(izz, i, 2) * sk(izz, i, 4);
            t(izz, 9) = sk(izz, i, 3) * sk(izz, i, 4);
            //
            // average point in k=constant face
            //
            t(izz, 4) = x(izz+1, 1, i) + x(izz+1, 1, i+1) + x(izz, 1, i) + x(izz, 1, i+1);
            t(izz, 5) = y(izz+1, 1, i) + y(izz+1, 1, i+1) + y(izz, 1, i) + y(izz, 1, i+1);
            t(izz, 6) = z(izz+1, 1, i) + z(izz+1, 1, i+1) + z(izz, 1, i) + z(izz, 1, i+1);
        }

        n = n - jdim;

        for (int izz = 1; izz <= n; izz++) {
            vol(izz, 1, i) = vol(izz, 1, i)
                           + (2.e0 * t(izz,       4) - t(izz, 10)) * t(izz,       7)
                           + (2.e0 * t(izz,       5) - t(izz, 11)) * t(izz,       8)
                           + (2.e0 * t(izz,       6) - t(izz, 12)) * t(izz,       9)
                           - (2.e0 * t(izz+jdim,  4) - t(izz, 10)) * t(izz+jdim,  7)
                           - (2.e0 * t(izz+jdim,  5) - t(izz, 11)) * t(izz+jdim,  8)
                           - (2.e0 * t(izz+jdim,  6) - t(izz, 12)) * t(izz+jdim,  9);
        }

        //
        // directed area and magnitude (i=constant face)
        //
        l = i + 1;
        n = jdim * kdim1 - 1;

        for (int izz = 1; izz <= n; izz++) {
            t(izz, 7) = si(izz, l, 1) * si(izz, l, 4);
            t(izz, 8) = si(izz, l, 2) * si(izz, l, 4);
            t(izz, 9) = si(izz, l, 3) * si(izz, l, 4);
            //
            // average point in i=constant face
            //
            t(izz, 13) = t(izz, 10) - t(izz, 13);
            t(izz, 14) = t(izz, 11) - t(izz, 14);
            t(izz, 15) = t(izz, 12) - t(izz, 15);

            vol(izz, 1, i) = vol(izz, 1, i)
                           - (2.e0 * t(izz, 13) - t(izz, 10)) * t(izz, 7)
                           - (2.e0 * t(izz, 14) - t(izz, 11)) * t(izz, 8)
                           - (2.e0 * t(izz, 15) - t(izz, 12)) * t(izz, 9);
            vol(izz, 1, i) = t1 * vol(izz, 1, i);
        }

        for (int kk = 1; kk <= kdim1; kk++) {
            vol(jdim, kk, i) = vol(jdim1, kk, i);
        }

        for (int izz = 1; izz <= jdim; izz++) {
            vol(izz, kdim, i) = vol(izz, kdim1, i);
        }
    } // end do 1000

    //
    // Negative volume check
    //
    if (negvol == 0) {
        for (i = 1; i <= idim1; i++) {
            for (int k = 1; k <= kdim; k++) {
                for (int j = 1; j <= jdim; j++) {
                    if ((float)vol(j, k, i) > 0.0f) continue;
                    nou(1) = std::min(nou(1) + 1, ibufdim);
                    std::snprintf(bou(nou(1), 1), 120,
                        "stopping ... negative volume at i,j,k=%5d%5d%5d block%5d",
                        i, j, k, nbl);
                    int m1 = -1;
                    termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                }
            }
        }
    } else {
        imin = 100000000;
        imax = 0;
        jmin = 100000000;
        jmax = 0;
        kmin = 100000000;
        kmax = 0;
        for (i = 1; i <= idim1; i++) {
            for (int k = 1; k <= kdim; k++) {
                for (int j = 1; j <= jdim; j++) {
                    if ((float)vol(j, k, i) < (float)eps) {
                        if (i < imin) imin = i;
                        if (i > imax) imax = i;
                        if (j < jmin) jmin = j;
                        if (j > jmax) jmax = j;
                        if (k < kmin) kmin = k;
                        if (k > kmax) kmax = k;
                        iflagv = 1;
                        vol(j, k, i) = ccomplex_ns::ccabs(vol(j, k, i));
                        if (iflagv1 == 1) {
                            nou(1) = std::min(nou(1) + 1, ibufdim);
                            std::snprintf(bou(nou(1), 1), 120,
                                "WARNING  ... negative volume at i,j,k=%5d%5d%5d block%5d not stopping!",
                                i, j, k, nbl);
                        }
                    }
                }
            }
        }
        imin = imin - 10;
        jmin = jmin - 10;
        kmin = kmin - 10;
        imax = imax + 10;
        jmax = jmax + 10;
        kmax = kmax + 10;
        if (imin < 1)    imin = 1;
        if (jmin < 1)    jmin = 1;
        if (kmin < 1)    kmin = 1;
        if (imax > idim) imax = idim;
        if (jmax > jdim) jmax = jdim;
        if (kmax > kdim) kmax = kdim;
    }

    return;
}

} // namespace cellvol_ns
