// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
// CFL3D - collx module
// Purpose: Restrict x, y, and z values to coarser meshes.

#include "collx.h"

namespace collx_ns {

void collx(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
           FortranArray3DRef<double> xx, FortranArray3DRef<double> yy, FortranArray3DRef<double> zz,
           int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2)
{
    int ii, iinc, i, kk, k, jj, j;

    ii   = 0;
    iinc = 2;
    if (idim == 2) iinc = 1;
    for (i = 1; i <= idim; i += iinc) {
        ii = ii + 1;
        kk = 0;
        for (k = 1; k <= kdim; k += 2) {
            kk = kk + 1;
            jj = 0;
            for (j = 1; j <= jdim; j += 2) {
                jj = jj + 1;
                xx(jj, kk, ii) = x(j, k, i);
                yy(jj, kk, ii) = y(j, k, i);
                zz(jj, kk, ii) = z(j, k, i);
            }
        }
    }
    return;
}

} // namespace collx_ns
