// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include <cstdio>
#include <cmath>
#include "l2norm.h"
#include "q8sdot.h"
#include <algorithm>

namespace l2norm_ns {

void l2norm(int& nbl, int& ntime, double& rmsl, int& irdq, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray3DRef<double> vol)
{
    int idim1, jdim1, kdim1, nplq, npl;
    int i, l, ipl, ii, j, k;
    int n, meq;

    idim1 = idim - 1;
    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    nplq  = std::min(idim1, 999000 / (jdim * kdim));
    npl   = nplq;
    rmsl  = 0.0;
    if (irdq == 1) goto label1000;

    for (i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        for (l = 1; l <= 5; l++) {
            for (ipl = 1; ipl <= npl; ipl++) {
                ii = i + ipl - 1;
                for (j = 1; j <= jdim; j++) {
                    res(j, kdim, ii, l) = 0.0;
                }
                for (k = 1; k <= kdim1; k++) {
                    res(jdim, k, ii, l) = 0.0;
                }
            }
        }
    }

    meq = 1;
    npl = nplq;
    for (i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        n = npl * jdim * kdim - jdim - 1;
        for (l = 1; l <= meq; l++) {
            FortranArray1DRef<double> res_slice(&res(1, 1, i, l), n);
            rmsl = rmsl + q8sdot_ns::q8sdot(n, res_slice, n, res_slice);
        }
    }

    return;

label1000:
    n = jdim * kdim - jdim - 1;
    for (i = 1; i <= idim1; i++) {
        rmsl = 0.0;
        for (l = 1; l <= 5; l++) {
            for (k = 1; k <= kdim1; k++) {
                res(jdim, k, i, l) = 0.0;
            }
            FortranArray1DRef<double> res_slice(&res(1, 1, i, l), n);
            rmsl = rmsl + q8sdot_ns::q8sdot(n, res_slice, n, res_slice);
        }
    }
    return;
}

void l2normall(int& nbl, int& ntime, FortranArray1DRef<double> rmsl, int& irdq, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray3DRef<double> vol)
{
    int idim1, jdim1, kdim1, nplq, npl;
    int i, l, ipl, ii, j, k;
    int n, meq;

    idim1 = idim - 1;
    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    nplq  = std::min(idim1, 999000 / (jdim * kdim));
    npl   = nplq;
    // rmsl = 0.e0 (array assignment)
    for (l = 1; l <= 5; l++) rmsl(l) = 0.0;
    if (irdq == 1) goto label1000;

    for (i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        for (l = 1; l <= 5; l++) {
            for (ipl = 1; ipl <= npl; ipl++) {
                ii = i + ipl - 1;
                for (j = 1; j <= jdim; j++) {
                    res(j, kdim, ii, l) = 0.0;
                }
                for (k = 1; k <= kdim1; k++) {
                    res(jdim, k, ii, l) = 0.0;
                }
            }
        }
    }

    meq = 5;
    npl = nplq;
    for (i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        n = npl * jdim * kdim - jdim - 1;
        for (l = 1; l <= meq; l++) {
            FortranArray1DRef<double> res_slice(&res(1, 1, i, l), n);
            rmsl(l) = rmsl(l) + q8sdot_ns::q8sdot(n, res_slice, n, res_slice);
        }
    }

    return;

label1000:
    n = jdim * kdim - jdim - 1;
    for (i = 1; i <= idim1; i++) {
        // rmsl = 0.e0 (array assignment inside loop)
        for (l = 1; l <= 5; l++) rmsl(l) = 0.0;
        for (l = 1; l <= 5; l++) {
            for (k = 1; k <= kdim1; k++) {
                res(jdim, k, i, l) = 0.0;
            }
            FortranArray1DRef<double> res_slice(&res(1, 1, i, l), n);
            rmsl(l) = rmsl(l) + q8sdot_ns::q8sdot(n, res_slice, n, res_slice);
        }
    }
    return;
}

} // namespace l2norm_ns
