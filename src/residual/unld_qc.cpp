// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
// CFL3D - unld_qc module
// Translated from Fortran to C++
//
// Purpose: Install the qtemp array used for message passing into
// the appropriate section of the qc array. Only the cell-center
// locations are unloaded, not the fill in locations on the idim,
// jdim and kdim planes.

#include "unld_qc.h"

namespace unld_qc_ns {

void unld_qc(FortranArray4DRef<double> qc, FortranArray4DRef<double> qtemp,
             int& jc, int& kc, int& ic, int& is, int& ie,
             int& js, int& je, int& ks, int& ke, int& ldim)
{
    int l, j, k, i;
    int jj, kk, ii;

    for (l = 1; l <= ldim; l++) {
        jj = js - 1;
        for (j = 1; j <= je - js; j++) {
            jj = jj + 1;
            kk = ks - 1;
            for (k = 1; k <= ke - ks; k++) {
                kk = kk + 1;
                ii = is - 1;
                for (i = 1; i <= ie - is; i++) {
                    ii = ii + 1;
                    qc(jj, kk, ii, l) = qtemp(j, k, i, l);
                }
            }
        }
    }
}

} // namespace unld_qc_ns
