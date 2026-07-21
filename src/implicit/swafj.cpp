// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "swafj.h"
#include "bsub.h"
#include "bsubp.h"

namespace swafj_ns {

void swafj(int& i, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> aj, FortranArray4DRef<double> bj,
           FortranArray4DRef<double> cj, FortranArray3DRef<double> f,
           int& nvt, FortranArray4DRef<double> res, int& iperd,
           FortranArray4DRef<double> gj, FortranArray4DRef<double> hj)
{
    // j-implicit k-sweep line inversions af

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;

    // load rhs (-residual) into f

    int kv = npl * kdim1;
    for (int ipl = 1; ipl <= npl; ipl++) {
        int ii  = i + ipl - 1;
        int kv0 = (ipl - 1) * kdim1;
        for (int l = 1; l <= 5; l++) {
            for (int k = 1; k <= kdim1; k++) {
                int k0 = k + kv0;

                int jj = 1 - kv;
                for (int jjj = 1; jjj <= jdim1; jjj++) {
                    jj = jj + kv;
                    f(k0 + jj - 1, 1, l) = res(jjj, k, ii, l);
                }
            }
        }
    }

    int n = kv * jdim1;
    for (int l = 1; l <= 5; l++) {
        for (int izz = 1; izz <= n; izz++) {
            f(izz, 1, l) = -f(izz, 1, l);
        }
    }

    // solve matrix equation

    int il  = 1;
    int iu  = jdim1;
    n       = kv;

    int id1 = npl * (kdim - 1);
    int i1  = 1;
    int i2  = n;
    if (iperd == 1) {
        bsubp_ns::bsubp(id1, jdim, aj, bj, cj, f, i1, i2, il, iu, gj, hj);
    } else {
        bsub_ns::bsub(id1, jdim, aj, bj, cj, f, i1, i2, il, iu);
    }

    // update delta q

    for (int ipl = 1; ipl <= npl; ipl++) {
        int ii  = i + ipl - 1;
        int kv0 = (ipl - 1) * kdim1;
        for (int l = 1; l <= 5; l++) {
            for (int k = 1; k <= kdim1; k++) {
                int k0 = k + kv0;

                int jj = 1 - kv;
                for (int jjj = 1; jjj <= jdim1; jjj++) {
                    jj = jj + kv;
                    res(jjj, k, ii, l) = f(k0 + jj - 1, 1, l);
                }
            }
        }
    }
    // label 2000 continue (no-op)
}

} // namespace swafj_ns
