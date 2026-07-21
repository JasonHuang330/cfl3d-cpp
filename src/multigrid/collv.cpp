// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
// CFL3D - collv module
// Purpose: Restrict volumes to coarser meshes.

#include "collv.h"

namespace collv_ns {

void collv(FortranArray3DRef<double> vol, FortranArray3DRef<double> volc,
           int& nj, int& nk, int& ni, int& jj2, int& kk2, int& ii2)
{
    int nj1, nk1, ni1, jjl, kkl, iil;
    int i, j, k, ii, jj, kk;

    nj1 = nj - 1;
    nk1 = nk - 1;
    ni1 = ni - 1;
    jjl = jj2 - 1;
    kkl = kk2 - 1;
    iil = ii2 - 1;
    ii  = 0;

    if (ni > 2) {
        for (i = 1; i <= ni1; i += 2) {
            ii = ii + 1;
            kk = 0;
            for (k = 1; k <= nk1; k += 2) {
                kk = kk + 1;
                jj = 0;
                for (j = 1; j <= nj1; j += 2) {
                    jj = jj + 1;
                    volc(jj, kk, ii) = vol(j,   k,   i  ) + vol(j+1, k,   i  )
                                     + vol(j,   k+1, i  ) + vol(j+1, k+1, i  )
                                     + vol(j,   k,   i+1) + vol(j+1, k,   i+1)
                                     + vol(j,   k+1, i+1) + vol(j+1, k+1, i+1);
                }
            }
        }
    } else {
        i  = 1;
        ii = 1;
        kk = 0;
        for (k = 1; k <= nk1; k += 2) {
            kk = kk + 1;
            jj = 0;
            for (j = 1; j <= nj1; j += 2) {
                jj = jj + 1;
                volc(jj, kk, ii) = vol(j,   k,   i) + vol(j+1, k,   i)
                                 + vol(j,   k+1, i) + vol(j+1, k+1, i);
            }
        }
    }

    for (k = 1; k <= kkl; k++) {
        for (i = 1; i <= iil; i++) {
            volc(jj2, k, i) = volc(jjl, k, i);
        }
    }

    for (i = 1; i <= iil; i++) {
        for (j = 1; j <= jj2; j++) {
            volc(j, kk2, i) = volc(j, kkl, i);
        }
    }

    return;
}

} // namespace collv_ns
