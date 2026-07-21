// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ld_qc.h"

namespace ld_qc_ns {

void ld_qc(FortranArray4DRef<double> qc, FortranArray4DRef<double> qtemp, int& jc, int& kc, int& ic, int& is, int& ie, int& js, int& je, int& ks, int& ke, int& ldim)
{
    // Purpose: Install the embedded subset of the qc array into qtemp
    // for message passing.

    int l, j, k, i, jj, kk, ii;

    for (l = 1; l <= ldim; l++) {
        jj = 0;
        for (j = js; j <= je; j++) {
            jj = jj + 1;
            kk = 0;
            for (k = ks; k <= ke; k++) {
                kk = kk + 1;
                ii = 0;
                for (i = is; i <= ie; i++) {
                    ii = ii + 1;
                    qtemp(jj, kk, ii, l) = qc(j, k, i, l);
                }
            }
        }
    }
}

} // namespace ld_qc_ns
