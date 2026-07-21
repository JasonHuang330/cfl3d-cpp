// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "hfluxl.h"
#include "dfluxpm.h"
#include "dfhat.h"

namespace hfluxl_ns {

void hfluxl(int& i, int& npl, double& xkap, int& idf, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qk0, FortranArray3DRef<double> sk,
            FortranArray3DRef<double> dhp, FortranArray3DRef<double> dhm,
            FortranArray2DRef<double> t, int& nvtq)
{
    // COMMON /fvfds/ rkap0(3), ifds(3)
    // ifds(3) -> cmn_fvfds.ifds[2]

    // h-flux-p-m jacobians

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int jv    = jdim1 * npl;
    int js    = jv + 1;

    // Loop 8008: fill q(j,kdim,ii,l) = qk0(j,ii,l,3)
    for (int ipl = 1; ipl <= npl; ipl++) {
        int ii = i + ipl - 1;
        for (int l = 1; l <= 5; l++) {
            for (int j = 1; j <= jdim; j++) {
                q(j, kdim, ii, l) = qk0(j, ii, l, 3);
            }
        }
    }

    // Loop 8000: fill t columns 36-39 (sk data) and t column 20 (sk(:,:,5))
    for (int k = 1; k <= kdim; k++) {
        int jkv = (k - 1) * jdim1 * npl + 1;
        int jk  = (k - 1) * jdim + 1;
        for (int ipl = 1; ipl <= npl; ipl++) {
            int jkv1 = jkv + (ipl - 1) * jdim1;
            int ii   = i + ipl - 1;
            for (int l = 1; l <= 4; l++) {
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + jkv1 - 1, 35 + l) = sk(izz + jk - 1, ii, l);
                }
            }
            for (int izz = 1; izz <= jdim1; izz++) {
                t(izz + jkv1 - 1, 20) = sk(izz + jk - 1, ii, 5);
            }
        }
    }

    int n = npl * jdim1 * kdim;

    // Loop 40: fill t columns 21-25 (q left) and 26-30 (q right)
    for (int k = 1; k <= kdim; k++) {
        int jkv = (k - 1) * jdim1 * npl + 1;
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii   = i + ipl - 1;
            int jkv1 = jkv + (ipl - 1) * jdim1;
            for (int l = 1; l <= 5; l++) {
                if (k > 1) {
                    for (int izz = 1; izz <= jdim1; izz++) {
                        t(izz + jkv1 - 1, 20 + l) = q(izz, k - 1, ii, l);
                    }
                } else {
                    for (int izz = 1; izz <= jdim1; izz++) {
                        t(izz + jkv1 - 1, 20 + l) = qk0(izz, ii, l, 1);
                    }
                }
                for (int izz = 1; izz <= jdim1; izz++) {
                    t(izz + jkv1 - 1, 25 + l) = q(izz, k, ii, l);
                }
            }
        }
    }

    // Build FortranArrayRef views into t columns for dfluxpm/dfhat calls
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

    // dhp: ipm = +1
    if (cmn_fvfds.ifds[2] == 0) {
        dfluxpm_ns::dfluxpm(t_col1, t_col2,
                            t_col36, t_col37, t_col38, t_col39,
                            t_col20, t_col21,
                            dhp, n, n, nvtq, ipm_pos);
    } else {
        dfhat_ns::dfhat(t_col36, t_col37, t_col38, t_col39,
                        t_col20, t_col21,
                        dhp, n, nvtq, ipm_pos);
    }

    // dhm: ipm = -1
    if (cmn_fvfds.ifds[2] == 0) {
        dfluxpm_ns::dfluxpm(t_col1, t_col2,
                            t_col36, t_col37, t_col38, t_col39,
                            t_col20, t_col26,
                            dhm, n, n, nvtq, ipm_neg);
    } else {
        dfhat_ns::dfhat(t_col36, t_col37, t_col38, t_col39,
                        t_col20, t_col26,
                        dhm, n, nvtq, ipm_neg);
    }
}

} // namespace hfluxl_ns
