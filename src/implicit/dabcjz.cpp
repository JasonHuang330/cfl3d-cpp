// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dabcjz.h"

namespace dabcjz_ns {

void dabcjz(int& i, int& npl, int& jdim, int& kdim, int& idim,
            FortranArray3DRef<double> a, FortranArray3DRef<double> b,
            FortranArray3DRef<double> c, FortranArray3DRef<double> blank)
{
    int j, ipl, k, ii;

    for (j = 1; j <= jdim; j++) {
        for (ipl = 1; ipl <= npl; ipl++) {
            ii = i + ipl - 1;
            for (k = 1; k <= kdim; k++) {
                a(k, ipl, j) =  a(k, ipl, j) * blank(j, k, ii);
                c(k, ipl, j) =  c(k, ipl, j) * blank(j, k, ii);
                b(k, ipl, j) = (b(k, ipl, j) * blank(j, k, ii)) + (1.0 - blank(j, k, ii));
            }
        }
    }
}

} // namespace dabcjz_ns
