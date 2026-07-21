// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "abckz.h"

namespace abckz_ns {

void abckz(int& i, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray5DRef<double> ak,
           FortranArray5DRef<double> bk,
           FortranArray5DRef<double> ck,
           FortranArray3DRef<double> blank)
{
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;

    // Zero out off-diagonal and diagonal elements of ak and ck
    // for blanked points (blank=0)
    for (int m = 1; m <= 5; m++) {
    for (int n = 1; n <= 5; n++) {
    for (int k = 1; k <= kdim1; k++) {
    for (int lpl = 1; lpl <= npl; lpl++) {
        int ii = i + lpl - 1;
        for (int j = 1; j <= jdim1; j++) {
            ak(j,lpl,k,m,n) = ak(j,lpl,k,m,n) * blank(j,k,ii);
            ck(j,lpl,k,m,n) = ck(j,lpl,k,m,n) * blank(j,k,ii);
        }
    }
    }
    }
    }

    // Modify diagonal elements of bk:
    // for blanked points (blank=0), set diagonal to 1.0
    for (int m = 1; m <= 5; m++) {
    for (int k = 1; k <= kdim1; k++) {
    for (int lpl = 1; lpl <= npl; lpl++) {
        int ii = i + lpl - 1;
        for (int j = 1; j <= jdim1; j++) {
            bk(j,lpl,k,m,m) = (bk(j,lpl,k,m,m) * blank(j,k,ii))
                             + (1.0 - blank(j,k,ii));
        }
    }
    }
    }
}

} // namespace abckz_ns
