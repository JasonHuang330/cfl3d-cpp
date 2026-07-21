// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "abciz.h"

namespace abciz_ns {

void abciz(int& k, int& npl, int& jdim, int& kdim, int& idim,
           FortranArray5DRef<double> ai, FortranArray5DRef<double> bi,
           FortranArray5DRef<double> ci, FortranArray3DRef<double> blank)
{
    int idim1 = idim - 1;
    int jdim1 = jdim - 1;

    // Zero off-diagonal elements (ai, ci) scaled by blank
    for (int m = 1; m <= 5; m++) {
        for (int n = 1; n <= 5; n++) {
            for (int i = 1; i <= idim1; i++) {
                for (int lpl = 1; lpl <= npl; lpl++) {
                    int kk = k + lpl - 1;
                    // cdir$ ivdep
                    for (int j = 1; j <= jdim1; j++) {
                        ai(j, lpl, i, m, n) = ai(j, lpl, i, m, n) * blank(j, kk, i);
                        ci(j, lpl, i, m, n) = ci(j, lpl, i, m, n) * blank(j, kk, i);
                    }
                }
            }
        }
    }

    // Modify diagonal elements of bi
    for (int m = 1; m <= 5; m++) {
        for (int i = 1; i <= idim1; i++) {
            for (int lpl = 1; lpl <= npl; lpl++) {
                int kk = k + lpl - 1;
                // cdir$ ivdep
                for (int j = 1; j <= jdim1; j++) {
                    bi(j, lpl, i, m, m) =   (bi(j, lpl, i, m, m) * blank(j, kk, i))
                                           + (1.0 - blank(j, kk, i));
                }
            }
        }
    }
}

} // namespace abciz_ns
