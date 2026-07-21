// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "colldat.h"

namespace colldat_ns {

void colldat(FortranArray4DRef<double> bcdata, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdatac, int& mdimc, int& ndimc)
{
    int l, m, n;
    int mm, nn, minc, ninc, mp, np, nnn, mmm;

    for (l = 1; l <= 12; l++) {
        mm   = 0;
        minc = 2;
        mp   = 1;
        if (mdim == 1) {
            minc = 1;
            mp   = 0;
        }
        for (m = 1; m <= mdim; m += minc) {
            mm = mm + 1;
            mm = std::min(mm, mdimc);
            nn   = 0;
            ninc = 2;
            np   = 1;
            if (ndim == 1) {
                ninc = 1;
                np   = 0;
            }
            for (n = 1; n <= ndim; n += ninc) {
                nn  = nn + 1;
                nn  = std::min(nn, ndimc);
                nnn = n + np;
                nnn = std::min(nnn, ndim);
                mmm = m + mp;
                mmm = std::min(mmm, mdim);
                bcdatac(mm, nn, 1, l) =
                    .25 * (bcdata(m, n, 1, l) + bcdata(mmm, n, 1, l) +
                           bcdata(m, nnn, 1, l) + bcdata(mmm, nnn, 1, l));
                bcdatac(mm, nn, 2, l) =
                    .25 * (bcdata(m, n, 2, l) + bcdata(mmm, n, 2, l) +
                           bcdata(m, nnn, 2, l) + bcdata(mmm, nnn, 2, l));
            }
        }
    }
}

} // namespace colldat_ns
