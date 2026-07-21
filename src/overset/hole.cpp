// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "hole.h"

namespace hole_ns {

void hole(int& i, int& npl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> res, FortranArray3DRef<double> blank)
{
//***********************************************************************
//     Purpose:  Zero out the rhs residuals for the blanked points.
//***********************************************************************

    int nn = jdim * kdim;
    for (int ipl = 1; ipl <= npl; ipl++) {
        int ii = i + ipl - 1;
        for (int l = 1; l <= 5; l++) {
            for (int n = 1; n <= nn; n++) {
                res(n, 1, ii, l) = res(n, 1, ii, l) * blank(n, 1, ii);
            }
        }
    }
}

} // namespace hole_ns
