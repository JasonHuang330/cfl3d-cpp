// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "abcjz.h"

namespace abcjz_ns {

void abcjz(int& i, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray5DRef<double> aj, FortranArray5DRef<double> bj,
           FortranArray5DRef<double> cj, FortranArray3DRef<double> blank)
{
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;

    for (int m = 1; m <= 5; m++) {
    for (int n = 1; n <= 5; n++) {
    for (int j = 1; j <= jdim1; j++) {
    for (int lpl = 1; lpl <= npl; lpl++) {
        int ii = i + lpl - 1;
        for (int k = 1; k <= kdim1; k++) {
            aj(k, lpl, j, m, n) = aj(k, lpl, j, m, n) * blank(j, k, ii);
            cj(k, lpl, j, m, n) = cj(k, lpl, j, m, n) * blank(j, k, ii);
        }
    }
    }
    }
    }

    for (int m = 1; m <= 5; m++) {
    for (int j = 1; j <= jdim1; j++) {
    for (int lpl = 1; lpl <= npl; lpl++) {
        int ii = i + lpl - 1;
        for (int k = 1; k <= kdim1; k++) {
            bj(k, lpl, j, m, m) = (bj(k, lpl, j, m, m) * blank(j, k, ii))
                                 + (1.0 - blank(j, k, ii));
        }
    }
    }
    }
}

} // namespace abcjz_ns
