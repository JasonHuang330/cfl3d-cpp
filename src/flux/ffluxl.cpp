// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ffluxl.h"
#include "dfluxpm.h"
#include "dfhat.h"

namespace ffluxl_ns {

void ffluxl(int& k, int& npl, double& xkap, int& idf, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qi0, FortranArray3DRef<double> si,
            FortranArray3DRef<double> dfp, FortranArray3DRef<double> dfm,
            FortranArray2DRef<double> t, int& nvtq)
{
    // f-flux-p-m jacobians

    int jdim1 = jdim - 1;
    int idim1 = idim - 1;
    int jv    = npl * jdim1;
    int js    = jv + 1;

    // Loop 8008: copy qi0(j,kk,l,3) -> q(j,kk,idim,l)
    for (int kpl = 1; kpl <= npl; kpl++) {
        int kk = k + kpl - 1;
        for (int l = 1; l <= 5; l++) {
            for (int j = 1; j <= jdim; j++) {
                q(j, kk, idim, l) = qi0(j, kk, l, 3);
            }
        }
    }

    // Loop 8000: fill t columns 36-39 (si l=1..4) and t col 6 (si l=5)
    for (int i = 1; i <= idim; i++) {
        int jiv = (i - 1) * jdim1 * npl + 1;
        for (int kpl = 1; kpl <= npl; kpl++) {
            int jiv1 = jiv + (kpl - 1) * jdim1;
            int kk   = k + kpl - 1;
            int ji   = (kk - 1) * jdim + 1;
            for (int l = 1; l <= 4; l++) {
                // cdir$ ivdep
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + jiv1 - 1, 35 + l) = si(izz + ji - 1, i, l);
                }
            }
            // cdir$ ivdep
            for (int izz = 1; izz <= jdim1; izz++) {
                t(izz + jiv1 - 1, 6) = si(izz + ji - 1, i, 5);
            }
        }
    }

    // create transposed q
    int n = npl * jdim1 * idim;
    for (int kpl = 1; kpl <= npl; kpl++) {
        int kk  = k + kpl - 1;
        int jiv = (kpl - 1) * jdim1 + 1;
        for (int l = 1; l <= 5; l++) {
            // cdir$ ivdep
            for (int izz = 1; izz <= jdim1; izz++) {
                t(izz + jiv - 1, 15 + l) = qi0(izz, kk, l, 1);
            }
            for (int i = 1; i <= idim; i++) {
                int jiv1 = (i - 1) * jdim1 * npl + jiv;
                // cdir$ ivdep
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + jiv1 - 1, 25 + l) = q(izz, kk, i, l);
                }
            }
        }
    }

    int n0 = n - jv;

    // first order
    for (int l = 1; l <= 5; l++) {
        // cdir$ ivdep
        for (int izz = 1; izz <= n0; izz++) {
            t(izz + js - 1, 20 + l) = t(izz, 25 + l);
        }
        // cdir$ ivdep
        for (int izz = 1; izz <= jv; izz++) {
            t(izz, 20 + l) = t(izz, 15 + l);
        }
    }

    // put unsteady metrics into t(20)
    // cdir$ ivdep
    for (int izz = 1; izz <= n; izz++) {
        t(izz, 20) = t(izz, 6);
    }

    // Build column slice refs for dfluxpm/dfhat calls
    // t(1,1): s1 (1D, nvtq elements)
    FortranArray1DRef<double> t_col1(&t(1, 1),   nvtq, 1);
    // t(1,2): s2 (1D, nvtq elements)
    FortranArray1DRef<double> t_col2(&t(1, 2),   nvtq, 1);
    // t(1,36): ax
    FortranArray1DRef<double> t_col36(&t(1, 36), nvtq, 1);
    // t(1,37): ay
    FortranArray1DRef<double> t_col37(&t(1, 37), nvtq, 1);
    // t(1,38): az
    FortranArray1DRef<double> t_col38(&t(1, 38), nvtq, 1);
    // t(1,39): area
    FortranArray1DRef<double> t_col39(&t(1, 39), nvtq, 1);
    // t(1,20): at
    FortranArray1DRef<double> t_col20(&t(1, 20), nvtq, 1);
    // t(1,21): q for dfp (2D, nvtq x 5)
    FortranArray2DRef<double> t_col21(&t(1, 21), nvtq, 5);
    // t(1,26): q for dfm (2D, nvtq x 5)
    FortranArray2DRef<double> t_col26(&t(1, 26), nvtq, 5);

    int ipm_pos = +1;
    int ipm_neg = -1;

    if (cmn_fvfds.ifds[0] == 0) {
        dfluxpm_ns::dfluxpm(t_col1, t_col2,
                            t_col36, t_col37, t_col38, t_col39,
                            t_col20, t_col21,
                            dfp, n, n, nvtq, ipm_pos);
    } else {
        dfhat_ns::dfhat(t_col36, t_col37, t_col38, t_col39,
                        t_col20, t_col21,
                        dfp, n, nvtq, ipm_pos);
    }

    if (cmn_fvfds.ifds[0] == 0) {
        dfluxpm_ns::dfluxpm(t_col1, t_col2,
                            t_col36, t_col37, t_col38, t_col39,
                            t_col20, t_col26,
                            dfm, n, n, nvtq, ipm_neg);
    } else {
        dfhat_ns::dfhat(t_col36, t_col37, t_col38, t_col39,
                        t_col20, t_col26,
                        dfm, n, nvtq, ipm_neg);
    }
}

} // namespace ffluxl_ns
