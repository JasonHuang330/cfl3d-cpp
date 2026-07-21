// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dabckz.h"

namespace dabckz_ns {

void dabckz(int& i, int& npl, int& jdim, int& kdim, int& idim,
            FortranArray3DRef<double> a, FortranArray3DRef<double> b,
            FortranArray3DRef<double> c, FortranArray3DRef<double> blank)
{
    int jdim1 = jdim - 1;

    for (int k = 1; k <= kdim; k++) {
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii = i + ipl - 1;
            for (int j = 1; j <= jdim1; j++) {
                a(j, ipl, k) = a(j, ipl, k) * blank(j, k, ii);
                c(j, ipl, k) = c(j, ipl, k) * blank(j, k, ii);
                b(j, ipl, k) = (b(j, ipl, k) * blank(j, k, ii)) + (1.0 - blank(j, k, ii));
            }
        }
    }
}

} // namespace dabckz_ns
