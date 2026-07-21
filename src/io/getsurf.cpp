// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "getsurf.h"

namespace getsurf_ns {

void getsurf(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
             FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti,
             FortranArray3DRef<double> xnm1, FortranArray3DRef<double> ynm1, FortranArray3DRef<double> znm1,
             FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
             FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
             FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
             FortranArray1DRef<double> wk, int& nwork, int& nbl,
             int& idim, int& jdim, int& kdim,
             int& nsurf, int& nsurfb,
             FortranArray1DRef<int> nsegdfrm, int& maxbl,
             FortranArray2DRef<int> idfrmseg, int& maxsegdg)
{
    int ll, is, i, j, k, ii, jj, kk;

    ll = nsurf * 9;

    for (is = 1; is <= nsegdfrm(nbl); is++) {

        //
        //        put i-surface points in list
        //
        if (icsi(nbl, is) == icsf(nbl, is)) {

            i  = icsi(nbl, is);
            ii = 1;
            if (icsi(nbl, is) == idim) ii = 2;
            for (j = jcsi(nbl, is); j <= jcsf(nbl, is); j++) {
                for (k = kcsi(nbl, is); k <= kcsf(nbl, is); k++) {
                    wk(ll + 1) = x(j, k, i);
                    wk(ll + 2) = y(j, k, i);
                    wk(ll + 3) = z(j, k, i);
                    wk(ll + 4) = delti(j, k, 1, ii);
                    wk(ll + 5) = delti(j, k, 2, ii);
                    wk(ll + 6) = delti(j, k, 3, ii);
                    wk(ll + 7) = xnm1(j, k, i);
                    wk(ll + 8) = ynm1(j, k, i);
                    wk(ll + 9) = znm1(j, k, i);
                    ll         = ll + 9;
                    nsurfb     = nsurfb + 1;
                }
            }

        }

        //
        //        put j-surface points in list
        //
        if (jcsi(nbl, is) == jcsf(nbl, is)) {

            j  = jcsi(nbl, is);
            jj = 1;
            if (jcsi(nbl, is) == jdim) jj = 2;
            for (i = icsi(nbl, is); i <= icsf(nbl, is); i++) {
                for (k = kcsi(nbl, is); k <= kcsf(nbl, is); k++) {
                    wk(ll + 1) = x(j, k, i);
                    wk(ll + 2) = y(j, k, i);
                    wk(ll + 3) = z(j, k, i);
                    wk(ll + 4) = deltj(k, i, 1, jj);
                    wk(ll + 5) = deltj(k, i, 2, jj);
                    wk(ll + 6) = deltj(k, i, 3, jj);
                    wk(ll + 7) = xnm1(j, k, i);
                    wk(ll + 8) = ynm1(j, k, i);
                    wk(ll + 9) = znm1(j, k, i);
                    ll         = ll + 9;
                    nsurfb     = nsurfb + 1;
                }
            }

        }

        //
        //        put k-surface points in list
        //
        if (kcsi(nbl, is) == kcsf(nbl, is)) {

            k  = kcsi(nbl, is);
            kk = 1;
            if (kcsi(nbl, is) == kdim) kk = 2;
            for (j = jcsi(nbl, is); j <= jcsf(nbl, is); j++) {
                for (i = icsi(nbl, is); i <= icsf(nbl, is); i++) {
                    wk(ll + 1) = x(j, k, i);
                    wk(ll + 2) = y(j, k, i);
                    wk(ll + 3) = z(j, k, i);
                    wk(ll + 4) = deltk(j, i, 1, kk);
                    wk(ll + 5) = deltk(j, i, 2, kk);
                    wk(ll + 6) = deltk(j, i, 3, kk);
                    wk(ll + 7) = xnm1(j, k, i);
                    wk(ll + 8) = ynm1(j, k, i);
                    wk(ll + 9) = znm1(j, k, i);
                    ll         = ll + 9;
                    nsurfb     = nsurfb + 1;
                }
            }

        }

    }

    return;
}

} // namespace getsurf_ns
