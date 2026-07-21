// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dabciz.h"

namespace dabciz_ns {

void dabciz(int& k, int& npl, int& jdim, int& kdim, int& idim,
            FortranArray3DRef<double> a, FortranArray3DRef<double> b,
            FortranArray3DRef<double> c, FortranArray3DRef<double> blank)
{
    int jdim1 = jdim - 1;

    for (int i = 1; i <= idim; i++) {
        for (int ipl = 1; ipl <= npl; ipl++) {
            int kk = k + ipl - 1;
            for (int j = 1; j <= jdim1; j++) {
                a(j, ipl, i) =  a(j, ipl, i) * blank(j, kk, i);
                c(j, ipl, i) =  c(j, ipl, i) * blank(j, kk, i);
                b(j, ipl, i) = (b(j, ipl, i) * blank(j, kk, i)) + (1.0 - blank(j, kk, i));
            }
        }
    }
}

} // namespace dabciz_ns
