// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "gfluxl.h"
#include "dfluxpm.h"
#include "dfhat.h"

namespace gfluxl_ns {

void gfluxl(int& i, int& npl, double& xkap, int& idf, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray3DRef<double> sj,
            FortranArray3DRef<double> dgp, FortranArray3DRef<double> dgm,
            FortranArray2DRef<double> t, int& nvtq)
{
    // COMMON /fvfds/ rkap0(3), ifds(3)
    // float rkap0[3], int32_t ifds[3]  (0-based)

    // g-flux-p-m jacobians

    int kdim1 = kdim - 1;

    int n  = jdim * kdim1;
    int n0 = n - 1;

    // Loop 8008: fill q(jdim,k,ii,l) = qj0(k,ii,l,3)
    for (int ipl = 1; ipl <= npl; ipl++) {
        int ii = i + ipl - 1;
        for (int l = 1; l <= 5; l++) {
            for (int k = 1; k <= kdim; k++) {
                q(jdim, k, ii, l) = qj0(k, ii, l, 3);
            }
        }
    }

    // Loop 8000: fill t columns 36-39 from sj, t col 20 from sj col 5
    for (int ipl = 1; ipl <= npl; ipl++) {
        int jkv1 = 1 + (ipl - 1) * jdim * kdim1;
        int ii   = i + ipl - 1;
        for (int l = 1; l <= 4; l++) {
            for (int izz = 1; izz <= n; izz++) {
                t(izz + jkv1 - 1, 35 + l) = sj(izz, ii, l);
            }
        }
        for (int izz = 1; izz <= n; izz++) {
            t(izz + jkv1 - 1, 20) = sj(izz, ii, 5);
        }
    }

    // Loop 40: fill t cols 21-25 and 26-30 from q, with boundary corrections from qj0
    for (int ipl = 1; ipl <= npl; ipl++) {
        int ii  = i + ipl - 1;
        int jkv = (ipl - 1) * jdim * kdim1;
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= n0; izz++) {
                t(izz + 1 + jkv, 20 + l) = q(izz, 1, ii, l);
            }
            for (int izz = 1; izz <= n; izz++) {
                t(izz + jkv, 25 + l) = q(izz, 1, ii, l);
            }
            for (int k = 1; k <= kdim1; k++) {
                int jk1         = jkv + (k - 1) * jdim + 1;
                t(jk1, 20 + l)  = qj0(k, ii, l, 1);
                int jk2         = jkv + k * jdim;
                t(jk2, 25 + l)  = q(jdim, k, ii, l);
            }
        }
    }

    n     = npl * n;
    int jkpro = jdim * kdim1 * npl;

    // Build column-slice refs for dfluxpm / dfhat calls
    // t(1,col) → pointer to start of column col in t(nvtq,39)
    FortranArray1DRef<double> t_col1 (&t(1,  1), nvtq, 1);
    FortranArray1DRef<double> t_col2 (&t(1,  2), nvtq, 1);
    FortranArray1DRef<double> t_col36(&t(1, 36), nvtq, 1);
    FortranArray1DRef<double> t_col37(&t(1, 37), nvtq, 1);
    FortranArray1DRef<double> t_col38(&t(1, 38), nvtq, 1);
    FortranArray1DRef<double> t_col39(&t(1, 39), nvtq, 1);
    FortranArray1DRef<double> t_col20(&t(1, 20), nvtq, 1);
    FortranArray2DRef<double> t_col21(&t(1, 21), nvtq, 5);
    FortranArray2DRef<double> t_col26(&t(1, 26), nvtq, 5);

    int ipm_pos = +1;
    int ipm_neg = -1;

    if (cmn_fvfds.ifds[1] == 0) {
        dfluxpm_ns::dfluxpm(t_col1, t_col2,
                            t_col36, t_col37, t_col38, t_col39,
                            t_col20, t_col21,
                            dgp, n, jkpro, nvtq, ipm_pos);
    } else {
        dfhat_ns::dfhat(t_col36, t_col37, t_col38, t_col39,
                        t_col20, t_col21,
                        dgp, n, nvtq, ipm_pos);
    }

    if (cmn_fvfds.ifds[1] == 0) {
        dfluxpm_ns::dfluxpm(t_col1, t_col2,
                            t_col36, t_col37, t_col38, t_col39,
                            t_col20, t_col26,
                            dgm, n, jkpro, nvtq, ipm_neg);
    } else {
        dfhat_ns::dfhat(t_col36, t_col37, t_col38, t_col39,
                        t_col20, t_col26,
                        dgm, n, nvtq, ipm_neg);
    }
}

} // namespace gfluxl_ns
